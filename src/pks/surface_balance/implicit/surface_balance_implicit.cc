/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon, Adam Atchley, Satish Karra

   DOCUMENT ME
   Surface Energy Balance for Snow Surface and Ground Surface
   Calculates Energy flux, rate or water, and water temperature
   entering through the surface skin.  Snow surface energy balance
   is calculated at equilibrium with ground/surface water and Air.

   ------------------------------------------------------------------------- */

#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"
#include "surface_balance_implicit.hh"

namespace Amanzi {
namespace SurfaceBalance {

SurfaceBalanceImplicit::SurfaceBalanceImplicit(
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           Teuchos::ParameterList& FElist,
           const Teuchos::RCP<TreeVector>& solution) :
    PKPhysicalBDFBase(plist, FElist, solution),
    PKDefaultBase(plist, FElist, solution),
    modify_predictor_advance_(false)
{
  // set up additional primary variables -- this is very hacky...
  // -- surface energy source
  Teuchos::ParameterList& esource_sublist =
      FElist.sublist("surface-conducted_energy_source");
  esource_sublist.set("evaluator name", "surface-conducted_energy_source");
  esource_sublist.set("field evaluator type", "primary variable");

  // -- surface mass source
  Teuchos::ParameterList& wsource_sublist =
      FElist.sublist("surface-mass_source");
  wsource_sublist.set("evaluator name", "surface-mass_source");
  wsource_sublist.set("field evaluator type", "primary variable");

  // -- subsurface mass source for VaporFlux at cell center
  Teuchos::ParameterList& w_v_source_sublist =
    FElist.sublist("mass_source");
  w_v_source_sublist.set("evaluator name", "mass_source");
  w_v_source_sublist.set("field evaluator type", "primary variable");

  // -- surface energy temperature
  Teuchos::ParameterList& wtemp_sublist =
      FElist.sublist("surface-mass_source_temperature");
  wtemp_sublist.set("evaluator name", "surface-mass_source_temperature");
  wtemp_sublist.set("field evaluator type", "primary variable");

  // Derivatives for PC
  eval_derivatives_ = plist_->get<bool>("evaluate source derivatives", true);

  // min wind speed
  min_wind_speed_ = plist_->get<double>("minimum wind speed [m/s]?", 1.0);
  wind_speed_ref_ht_ = plist_->get<double>("wind speed reference height [m]", 2.0);

  // implicit/explicit snow precip
  implicit_snow_ = plist_->get<bool>("implicit snow precipitation", false);

  // Reading in Longwave Radation
  longwave_input_ = plist_->get<bool>("Longwave Input", false);

  // snow thermal conductivity model parameters
  snow_lambda_fresh_ = plist_->get<double>("thermal conductivity of fresh snow [W/m K]", 0.029);
  snow_lambda_exponent_ = plist_->get<double>("thermal conductivity aging exponent for snow", 2.0);

  // transition snow depth
  snow_ground_trans_ = plist_->get<double>("snow-ground transitional depth", 0.02);
  min_snow_trans_ = plist_->get<double>("minimum snow transitional depth", 1.e-8);
  if (min_snow_trans_ < 0. || snow_ground_trans_ < min_snow_trans_) {
    Errors::Message message("Invalid parameters: snow-ground transitional depth or minimum snow transitional depth.");
    Exceptions::amanzi_throw(message);
  }

  // modify predictor by calling advance -- this is cheap and sets up BCs
  // correctly for subsurface's call to ModifyPredictorConsistentFaces()
  modify_predictor_advance_ = plist_->get<bool>("modify predictor by advancing", false);
}

// main methods
// -- Setup data.
void
SurfaceBalanceImplicit::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);
  subsurf_mesh_ = S->GetMesh(); // needed for VPL, which is treated as subsurface source

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)->
      SetComponent("cell", AmanziMesh::CELL, 1);
  

  // requirements: other primary variables
  Teuchos::RCP<FieldEvaluator> fm;
  S->RequireField("surface-conducted_energy_source", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface-conducted_energy_source");
  fm = S->GetFieldEvaluator("surface-conducted_energy_source");
  pvfe_esource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_esource_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField("surface-mass_source", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface-mass_source");
  fm = S->GetFieldEvaluator("surface-mass_source");
  pvfe_wsource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_wsource_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField("mass_source", name_)->SetMesh(subsurf_mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("mass_source");
  fm = S->GetFieldEvaluator("mass_source");
  pvfe_w_v_source_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_w_v_source_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField("surface-mass_source_temperature", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface-mass_source_temperature");
  fm = S->GetFieldEvaluator("surface-mass_source_temperature");
  pvfe_wtemp_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_wtemp_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  // requirements: source derivatives
  if (eval_derivatives_) {
    S->RequireField("dsurface-conducted_energy_source_dsurface-temperature", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  }

  // requirements: diagnostic variables
  S->RequireField("albedo",name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField("albedo",name_)->set_io_checkpoint(false);
  S->RequireField("evaporative_flux",name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField("evaporative_flux",name_)->set_io_checkpoint(false);
  S->RequireField("qE_latent_heat",name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField("qE_latent_heat",name_)->set_io_checkpoint(false);
  S->RequireField("qE_sensible_heat",name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField("qE_sensible_heat",name_)->set_io_checkpoint(false);
  S->RequireField("qE_lw_out",name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField("qE_lw_out",name_)->set_io_checkpoint(false);
  
  // requirements: independent variables (data from MET)
  S->RequireFieldEvaluator("incoming_shortwave_radiation");
  S->RequireField("incoming_shortwave_radiation")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  if (longwave_input_) {
     S->RequireFieldEvaluator("incoming_longwave_radiation");
     S->RequireField("incoming_longwave_radiation")->SetMesh(mesh_)
        ->AddComponent("cell", AmanziMesh::CELL, 1);
  }

  S->RequireFieldEvaluator("air_temperature");
  S->RequireField("air_temperature")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("relative_humidity");
  S->RequireField("relative_humidity")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("wind_speed");
  S->RequireField("wind_speed")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("precipitation_rain");
  S->RequireField("precipitation_rain")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("precipitation_snow");
  S->RequireField("precipitation_snow")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirements: stored secondary variables
  S->RequireField("snow_density", name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField("snow_age", name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField("snow_temperature", name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireField("stored_SWE", name_)->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface-temperature");
  S->RequireField("surface-temperature")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface-pressure");
  S->RequireField("surface-pressure")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("ponded_depth");
  S->RequireField("ponded_depth")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

 //   S->RequireFieldEvaluator("saturation_liquid");
   S->RequireField("saturation_liquid")->SetMesh(subsurf_mesh_)
       ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("unfrozen_fraction");
  S->RequireField("unfrozen_fraction")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface-porosity");
  S->RequireField("surface-porosity")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

}

// -- Initialize owned (dependent) variables.
void
SurfaceBalanceImplicit::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::initialize(S);

  SEBPhysics::SEB seb;

  // initialize snow density, age
  ASSERT(plist_->isSublist("initial condition"));
  Teuchos::ParameterList& ic_list = plist_->sublist("initial condition");
  if (!S->GetField("snow_density")->initialized()) {
    if (ic_list.isParameter("restart file")) {
      // initialize density, age from restart file
      S->GetField("snow_density", name_)->Initialize(ic_list);
      S->GetField("snow_density", name_)->set_initialized();
      S->GetField("snow_age", name_)->Initialize(ic_list);
      S->GetField("snow_age", name_)->set_initialized();
    } else {
      // initialize density to fresh powder, age to 0
      S->GetFieldData("snow_density",name_)->PutScalar(seb.params.density_freshsnow);
      S->GetField("snow_density", name_)->set_initialized();
      S->GetFieldData("snow_age",name_)->PutScalar(0.);
      S->GetField("snow_age", name_)->set_initialized();
    }
  }

  // initialize swe consistently with snow height and density
  Epetra_MultiVector& swe = *S->GetFieldData("stored_SWE",name_)->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_ht = *S->GetFieldData(key_)->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_dens = *S->GetFieldData("snow_density")->ViewComponent("cell",false);
  for (int c=0; c!=swe.MyLength(); ++c) {
    swe[0][c] = snow_ht[0][c] * snow_dens[0][c] / seb.in.vp_ground.density_w;
  }
  S->GetField("stored_SWE", name_)->set_initialized();

  // initialize snow temp
  S->GetFieldData("snow_temperature",name_)->PutScalar(0.);
  S->GetField("snow_temperature", name_)->set_initialized();

  // initialize sources, temps
  S->GetFieldData("surface-conducted_energy_source",name_)->PutScalar(0.);
  S->GetField("surface-conducted_energy_source",name_)->set_initialized();

  if (eval_derivatives_) {
    S->GetFieldData("dsurface-conducted_energy_source_dsurface-temperature",name_)->PutScalar(0.);
    S->GetField("dsurface-conducted_energy_source_dsurface-temperature",name_)->set_initialized();
  }

  S->GetFieldData("surface-mass_source",name_)->PutScalar(0.);
  S->GetField("surface-mass_source",name_)->set_initialized();

  S->GetFieldData("mass_source",name_)->PutScalar(0.);
  S->GetField("mass_source",name_)->set_initialized();

  S->GetFieldData("surface-mass_source_temperature",name_)->PutScalar(273.15);
  S->GetField("surface-mass_source_temperature",name_)->set_initialized();

  // initialize diagnostics
  S->GetField("albedo",name_)->set_initialized();
  S->GetField("evaporative_flux",name_)->set_initialized();
  S->GetField("qE_latent_heat",name_)->set_initialized();
  S->GetField("qE_sensible_heat",name_)->set_initialized();
  S->GetField("qE_lw_out",name_)->set_initialized();
}


// computes the non-linear functional g = g(t,u,udot)
void
SurfaceBalanceImplicit::Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                            Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  Teuchos::OSTab tab = vo_->getOSTab();
  double dt = t_new - t_old;
  double T_eps = 1.e-10;

  bool debug = false;
  Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
  int rank = mesh_->get_comm()->MyPID();

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << dt << std::endl;
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("new snow depth");
    vecs.push_back(u_new->Data().ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();
  }

  // pull residual vector
  Epetra_MultiVector& res = *g->Data()->ViewComponent("cell",false);

  // pull old snow data
  const Epetra_MultiVector& snow_depth_old = *u_old->Data()->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_age_old = *S_inter_->GetFieldData("snow_age")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_dens_old = *S_inter_->GetFieldData("snow_density")
      ->ViewComponent("cell",false);
    const Epetra_MultiVector& stored_SWE_old = *S_inter_->GetFieldData("stored_SWE")
      ->ViewComponent("cell",false);

  // pull current snow data
  const Epetra_MultiVector& snow_depth_new = *u_new->Data()->ViewComponent("cell",false);
  Epetra_MultiVector& snow_temp_new = *S_next_->GetFieldData("snow_temperature", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& snow_age_new = *S_next_->GetFieldData("snow_age", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& snow_dens_new = *S_next_->GetFieldData("snow_density", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& stored_SWE_new = *S_next_->GetFieldData("stored_SWE", name_)
      ->ViewComponent("cell",false);

  // pull diagnostics
  Epetra_MultiVector& albedo = *S_next_->GetFieldData("albedo",name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& evaporative_flux = *S_next_->GetFieldData("evaporative_flux",name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& qE_latent_heat = *S_next_->GetFieldData("qE_latent_heat",name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& qE_sensible_heat = *S_next_->GetFieldData("qE_sensible_heat",name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& qE_lw_out = *S_next_->GetFieldData("qE_lw_out",name_)
      ->ViewComponent("cell",false);

  // pull ATS data
  S_next_->GetFieldEvaluator("surface-temperature")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_temp =
      *S_next_->GetFieldData("surface-temperature")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("surface-pressure")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_pres =
      *S_next_->GetFieldData("surface-pressure")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& ponded_depth =
      *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell", false);

 //  S_next_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_next_.ptr(), name_);
   const Epetra_MultiVector& saturation_liquid =
     *S_next_->GetFieldData("saturation_liquid")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("unfrozen_fraction")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& unfrozen_fraction =
      *S_next_->GetFieldData("unfrozen_fraction")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("surface-porosity")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_porosity =
      *S_next_->GetFieldData("surface-porosity")->ViewComponent("cell", false);


  // pull Met data
  S_next_->GetFieldEvaluator("air_temperature")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& air_temp =
      *S_next_->GetFieldData("air_temperature")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("incoming_shortwave_radiation")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& incoming_shortwave =
      *S_next_->GetFieldData("incoming_shortwave_radiation")->ViewComponent("cell", false);

 Teuchos::RCP<const Epetra_MultiVector> incoming_longwave = Teuchos::null;
  if (longwave_input_) {
  S_next_->GetFieldEvaluator("incoming_longwave_radiation")->HasFieldChanged(S_next_.ptr(), name_);
       incoming_longwave =
      S_next_->GetFieldData("incoming_longwave_radiation")->ViewComponent("cell", false);
  }

  S_next_->GetFieldEvaluator("relative_humidity")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& relative_humidity =
      *S_next_->GetFieldData("relative_humidity")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("wind_speed")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& wind_speed =
      *S_next_->GetFieldData("wind_speed")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("precipitation_rain")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& precip_rain =
      *S_next_->GetFieldData("precipitation_rain")->ViewComponent("cell", false);

  // snow precip need not be updated each iteration
  if (implicit_snow_) {
    S_next_->GetFieldEvaluator("precipitation_snow")->HasFieldChanged(S_next_.ptr(), name_);
  } else {
    S_inter_->GetFieldEvaluator("precipitation_snow")->HasFieldChanged(S_inter_.ptr(), name_);
  }
  const Epetra_MultiVector& precip_snow = implicit_snow_ ?
    *S_next_->GetFieldData("precipitation_snow")->ViewComponent("cell", false) :
    *S_inter_->GetFieldData("precipitation_snow")->ViewComponent("cell", false);

  // pull additional primary variable data
  Epetra_MultiVector& surf_energy_flux =
      *S_next_->GetFieldData("surface-conducted_energy_source", name_)->ViewComponent("cell", false);
  Teuchos::RCP<Epetra_MultiVector> dsurf_energy_flux_dT;
  if (eval_derivatives_) {
    dsurf_energy_flux_dT = S_next_->GetFieldData("dsurface-conducted_energy_source_dsurface-temperature", name_)
      ->ViewComponent("cell", false);
  }

  Epetra_MultiVector& surf_water_flux =
      *S_next_->GetFieldData("surface-mass_source", name_)->ViewComponent("cell", false);

  Epetra_MultiVector& vapor_flux =
      *S_next_->GetFieldData("mass_source", name_)->ViewComponent("cell", false);
  vapor_flux.PutScalar(0.);

  Epetra_MultiVector& surf_water_flux_temp =
      *S_next_->GetFieldData("surface-mass_source_temperature", name_)->ViewComponent("cell", false);


  unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (unsigned int c=0; c!=ncells; ++c) { // START CELL LOOP  ##########################
    dcvo = Teuchos::null;
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) dcvo = db_->GetVerboseObject(c, rank);
    Teuchos::OSTab dctab = dcvo == Teuchos::null ? vo_->getOSTab() : dcvo->getOSTab();

    SEBPhysics::SEB seb;

    double snow_depth = snow_depth_old[0][c];
    double snow_dens = snow_dens_old[0][c];
    double swe = stored_SWE_old[0][c];
    double snow_age = snow_age_old[0][c];
    if (!implicit_snow_) {
      double Ps = std::max(precip_snow[0][c], 0.);
      swe += dt * Ps;
      snow_depth += dt * Ps * seb.in.vp_ground.density_w / seb.params.density_freshsnow;
      snow_age = swe > 0. ? (snow_age * stored_SWE_old[0][c]) / swe : 0.; // age of new snow is 0 here, as increment happens internally
      snow_dens = snow_depth > 0. ? swe * seb.in.vp_ground.density_w / snow_depth : seb.params.density_freshsnow;
    }

    if (snow_depth >= snow_ground_trans_ ||
        snow_depth < min_snow_trans_) {
      // Evaluate the model as usual
      // Initialize the SEB object
      seb.in.dt = dt;

      // -- ground properties
      seb.in.vp_ground.temp = surf_temp[0][c];
      seb.in.vp_ground.pressure = surf_pres[0][c];
      AmanziMesh::Entity_ID subsurf_f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
      AmanziMesh::Entity_ID_List cells;
      subsurf_mesh_->face_get_cells(subsurf_f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);
      seb.in.surf.saturation_liquid = saturation_liquid[0][cells[0]];

      // -- snow properties
      seb.in.snow_old.ht = snow_depth < min_snow_trans_ ? 0. : snow_depth;
      seb.in.snow_old.density = snow_dens;
      seb.in.snow_old.age = snow_age;
      seb.in.snow_old.SWE = swe;
      seb.out.snow_new = seb.in.snow_old;
      seb.in.vp_snow.temp = 273.15;

      // -- snow thermal conductivity parameters
      seb.params.snoK = snow_lambda_fresh_;
      seb.params.snoKX = snow_lambda_exponent_;

      // -- met data
      seb.params.Zr = wind_speed_ref_ht_;
      seb.in.met.Us = std::max(wind_speed[0][c], min_wind_speed_);
      seb.in.met.QswIn = incoming_shortwave[0][c];
      seb.in.met.Ps = implicit_snow_ ? std::max(precip_snow[0][c],0.) : 0.; // protect against wayward snow distribution models
      seb.in.met.Pr = precip_rain[0][c];
      seb.in.met.vp_air.temp = air_temp[0][c];
      seb.in.met.vp_air.relative_humidity = relative_humidity[0][c];
     
     if (longwave_input_) {
          seb.in.met.QlwIn = (*incoming_longwave)[0][c];
      }else{     
          seb.in.met.vp_air.UpdateVaporPressure();
          double e_air = std::pow(10*seb.in.met.vp_air.actual_vaporpressure, seb.in.met.vp_air.temp / 2016.);
          e_air = 1.08 * (1 - std::exp(-e_air));
          seb.in.met.QlwIn = e_air * seb.params.stephB * std::pow(seb.in.met.vp_air.temp,4); // Add if statement here ~ AA
      }

      // -- smoothed/interpolated surface properties
      SEBPhysics::SurfaceParams surf_pars;

      SEBPhysics::Partition al_part = SEBPhysics::Partitioner()
          .CalcPartition(seb.in.snow_old.ht, ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb.in.surf.albedo = al_part.Interpolate(
          SEBPhysics::CalcAlbedoSnow(seb.in.snow_old.density),
          surf_pars.a_water, surf_pars.a_ice, surf_pars.a_tundra);

      SEBPhysics::Partition other_part = SEBPhysics::Partitioner(0.02, 0.02)
          .CalcPartition(seb.in.snow_old.ht, ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb.in.surf.emissivity = other_part.Interpolate(surf_pars.e_snow,
              surf_pars.e_water, surf_pars.e_ice, surf_pars.e_tundra);
      seb.in.vp_ground.porosity = other_part.Interpolate(1., 1., 1., surf_porosity[0][c]);

      // -- roughness factor
      seb.in.surf.Zo = SEBPhysics::CalcRoughnessFactor(seb.in.met.vp_air.temp);

      // Run the model
      SEBPhysics::CalculateSurfaceBalance(seb, false, dcvo);

      // Evaluate the residual
      res[0][c] =  snow_depth_new[0][c] - seb.out.snow_new.ht;

      // Pull the output
      // -- fluxes
      surf_energy_flux[0][c] = 1.e-6 * seb.out.eb.fQc; // convert to MJ for ATS
      surf_water_flux[0][c] = seb.out.mb.MWg;
      surf_water_flux_temp[0][c] = seb.out.mb.MWg_temp;

      // -- vapor flux to cells
      //     surface vapor flux is treated as a volumetric source for the subsurface.
      // surface mass sources are in m^3 water / (m^2 s)
      // subsurface mass sources are in mol water / (m^3 s)
      vapor_flux[0][cells[0]] = seb.out.mb.MWg_subsurf
          * mesh_->cell_volume(c) * seb.in.vp_ground.density_w / 0.0180153
          / subsurf_mesh_->cell_volume(cells[0]);

      // -- snow properties
      snow_age_new[0][c] = seb.out.snow_new.age;
      snow_dens_new[0][c] = seb.out.snow_new.density;
      snow_temp_new[0][c] = seb.in.vp_snow.temp;
      stored_SWE_new[0][c] = seb.out.snow_new.SWE;

      // -- diagnostics
      albedo[0][c] = seb.in.surf.albedo;
      evaporative_flux[0][c] = seb.out.mb.Me;
      qE_latent_heat[0][c] = seb.out.eb.fQe;
      qE_sensible_heat[0][c] = seb.out.eb.fQh;
      qE_lw_out[0][c] = seb.out.eb.fQlwOut;

      if (eval_derivatives_) {
        // evaluate FD derivative of energy flux wrt surface temperature
        SEBPhysics::SEB seb2(seb);
        seb.params.Zr = wind_speed_ref_ht_;
        seb2.in.vp_ground.temp += T_eps;
        // for now ignore the effect on unfrozen fraction, and therefore on albedo and emissivity
        SEBPhysics::CalculateSurfaceBalance(seb2);
        (*dsurf_energy_flux_dT)[0][c] = 1.e-6 * (seb2.out.eb.fQc - seb.out.eb.fQc) / T_eps; // MJ
      }

    } else {
      // Evaluate the model twice -- once as bare ground, once with snow, using
      // an area-averaged subgrid model to smooth between the two end-members.
      // The area-weighting parameter is theta:
      double theta = snow_depth / snow_ground_trans_;

      // Evaluate the model as usual
      // Initialize the SEB object
      seb.in.dt = dt;

      // -- ground properties
      seb.in.vp_ground.temp = surf_temp[0][c];
      seb.in.vp_ground.pressure = surf_pres[0][c];
      AmanziMesh::Entity_ID subsurf_f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
      AmanziMesh::Entity_ID_List cells;
      subsurf_mesh_->face_get_cells(subsurf_f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);
      seb.in.surf.saturation_liquid = saturation_liquid[0][cells[0]];

      // -- snow properties
      seb.in.snow_old.ht = snow_ground_trans_;
      seb.in.snow_old.density = snow_dens;
      seb.in.snow_old.age = snow_age;
      seb.in.snow_old.SWE = swe / theta;
      seb.out.snow_new = seb.in.snow_old;
      seb.in.vp_snow.temp = 273.15;

      // -- snow thermal conductivity parameters
      seb.params.snoK = snow_lambda_fresh_;
      seb.params.snoKX = snow_lambda_exponent_;

      // -- met data
      seb.params.Zr = wind_speed_ref_ht_;
      seb.in.met.Us = std::max(wind_speed[0][c], min_wind_speed_);
      seb.in.met.QswIn = incoming_shortwave[0][c];
      seb.in.met.Ps = implicit_snow_ ? std::max(precip_snow[0][c],0.) : 0.;
      seb.in.met.Pr = precip_rain[0][c];
      seb.in.met.vp_air.temp = air_temp[0][c];
      seb.in.met.vp_air.relative_humidity = relative_humidity[0][c];

      if (longwave_input_) {
	seb.in.met.QlwIn = (*incoming_longwave)[0][c];
      } else {
	seb.in.met.vp_air.UpdateVaporPressure();
	double e_air = std::pow(10*seb.in.met.vp_air.actual_vaporpressure, seb.in.met.vp_air.temp / 2016.);
	e_air = 1.08 * (1 - std::exp(-e_air));
	seb.in.met.QlwIn = e_air * seb.params.stephB * std::pow(seb.in.met.vp_air.temp,4); // Add if statement here ~ AA
      }

      // -- smoothed/interpolated surface properties
      SEBPhysics::SurfaceParams surf_pars;

      SEBPhysics::Partition al_part = SEBPhysics::Partitioner()
          .CalcPartition(seb.in.snow_old.ht, ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb.in.surf.albedo = al_part.Interpolate(
          SEBPhysics::CalcAlbedoSnow(seb.in.snow_old.density),
          surf_pars.a_water, surf_pars.a_ice, surf_pars.a_tundra);

      SEBPhysics::Partition other_part = SEBPhysics::Partitioner(0.02, 0.02)
          .CalcPartition(seb.in.snow_old.ht, ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb.in.surf.emissivity = other_part.Interpolate(surf_pars.e_snow,
              surf_pars.e_water, surf_pars.e_ice, surf_pars.e_tundra);
      seb.in.vp_ground.porosity = other_part.Interpolate(1., 1., 1., surf_porosity[0][c]);

      // -- roughness factor
      seb.in.surf.Zo = SEBPhysics::CalcRoughnessFactor(seb.in.met.vp_air.temp);

      // Initialize a second model, with bare ground
      SEBPhysics::SEB seb_bare(seb);
      seb_bare.params.Zr = wind_speed_ref_ht_;

      // -- snow properties
      seb_bare.in.snow_old.ht = 0.;
      seb_bare.out.snow_new.ht = 0.;

      // -- smoothed interpolated surface properties
      al_part = SEBPhysics::Partitioner()
          .CalcPartition(seb_bare.in.snow_old.ht, ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb_bare.in.surf.albedo = al_part.Interpolate(
          SEBPhysics::CalcAlbedoSnow(seb_bare.in.snow_old.density),
          surf_pars.a_water, surf_pars.a_ice, surf_pars.a_tundra);

      other_part = SEBPhysics::Partitioner(0.02, 0.02)
          .CalcPartition(seb_bare.in.snow_old.ht, ponded_depth[0][c], unfrozen_fraction[0][c]);
      seb_bare.in.surf.emissivity = other_part.Interpolate(surf_pars.e_snow,
              surf_pars.e_water, surf_pars.e_ice, surf_pars.e_tundra);
      seb_bare.in.vp_ground.porosity = other_part.Interpolate(1., 1., 1., surf_porosity[0][c]);

      // Run the model for both snowy case and bare case
      SEBPhysics::CalculateSurfaceBalance(seb, false, dcvo);
      SEBPhysics::CalculateSurfaceBalance(seb_bare, false, dcvo);

      // Evaluate the residual
      double snow_depth_new_tmp = theta * seb.out.snow_new.ht + (1-theta) * seb_bare.out.snow_new.ht;
      res[0][c] = snow_depth_new[0][c] - snow_depth_new_tmp;

      // Pull the output
      // -- fluxes
      surf_energy_flux[0][c] = 1.e-6 * (theta * seb.out.eb.fQc + (1-theta) * seb_bare.out.eb.fQc); // MJ
      surf_water_flux[0][c] = theta * seb.out.mb.MWg + (1-theta) * seb_bare.out.mb.MWg;
      if (surf_water_flux[0][c] > 0.) {
        if (seb.out.mb.MWg > 0.) {
          if (seb_bare.out.mb.MWg > 0.) {
            // both are positive, sources from both, average temp
            surf_water_flux_temp[0][c] = (theta * seb.out.mb.MWg * seb.out.mb.MWg_temp
                                          + (1-theta) * seb_bare.out.mb.MWg * seb_bare.out.mb.MWg_temp) / surf_water_flux[0][c];
          } else {
            // snow portion is positive, bare is negative, so take only from snow
            surf_water_flux_temp[0][c] = seb.out.mb.MWg_temp;
          }
        } else {
          // snow portion is negative, take only from snow.  Note if both are negative, both temps are just the air temp, which is fine.
          surf_water_flux_temp[0][c] = seb_bare.out.mb.MWg_temp;
        }
      } else {
        surf_water_flux_temp[0][c] = seb.out.mb.MWg_temp;
      }
      // -- vapor flux to cells
      //     surface vapor flux is treated as a volumetric source for the subsurface.
      // surface mass sources are in m^3 water / (m^2 s)
      // subsurface mass sources are in mol water / (m^3 s)
      double mean_flux = theta * seb.out.mb.MWg_subsurf + (1-theta) * seb_bare.out.mb.MWg_subsurf;
      vapor_flux[0][cells[0]] = mean_flux
          * mesh_->cell_volume(c) * seb.in.vp_ground.density_w / 0.0180153
          / subsurf_mesh_->cell_volume(cells[0]);

      // -- snow properties: SWE averaged
      double total_swe = theta * seb.out.snow_new.ht * seb.out.snow_new.density / seb.in.vp_ground.density_w
          + (1-theta) * seb_bare.out.snow_new.ht * seb_bare.out.snow_new.density / seb_bare.in.vp_ground.density_w;
      if (snow_depth_new_tmp > 0.) {
        snow_age_new[0][c] = (theta * seb.out.snow_new.ht * seb.out.snow_new.age
                              + (1-theta) * seb_bare.out.snow_new.ht * seb_bare.out.snow_new.age) / snow_depth_new_tmp;
        snow_dens_new[0][c] = total_swe * seb.in.vp_ground.density_w / snow_depth_new_tmp;
        stored_SWE_new[0][c] = total_swe;
      } else {
        snow_age_new[0][c] = 0.;
        snow_dens_new[0][0] = seb.out.snow_new.density;
        stored_SWE_new[0][c] = total_swe;
      }
      snow_temp_new[0][c] = seb.in.vp_snow.temp;

      // -- diagnostics
      albedo[0][c] = theta * seb.in.surf.albedo + (1-theta) * seb_bare.in.surf.albedo;
      evaporative_flux[0][c] = theta * seb.out.mb.Me + (1-theta) * seb.out.mb.Me;
      qE_latent_heat[0][c] = theta * seb.out.eb.fQe + (1-theta) * seb.out.eb.fQe;
      qE_sensible_heat[0][c] = theta * seb.out.eb.fQh + (1-theta) * seb.out.eb.fQh;
      qE_lw_out[0][c] = theta * seb.out.eb.fQlwOut + (1-theta) * seb.out.eb.fQlwOut;

      // Evaluate derivatives, if requested
      if (eval_derivatives_) {
        // evaluate FD derivative of energy flux wrt surface temperature
        SEBPhysics::SEB seb2(seb);
        seb2.params.Zr = wind_speed_ref_ht_;
        SEBPhysics::SEB seb2_bare(seb_bare);
        seb2_bare.params.Zr = wind_speed_ref_ht_;
        seb2.in.vp_ground.temp += T_eps;
        seb2_bare.in.vp_ground.temp += T_eps;
        // for now ignore the effect on unfrozen fraction, and therefore on albedo and emissivity
        SEBPhysics::CalculateSurfaceBalance(seb2);
        SEBPhysics::CalculateSurfaceBalance(seb2_bare);
        double eflux2 = 1.e-6 * (theta * seb2.out.eb.fQc + (1-theta) * seb2_bare.out.eb.fQc);
        
        (*dsurf_energy_flux_dT)[0][c] = (eflux2 - surf_energy_flux[0][c]) / T_eps; // MJ
      }
    }
  }  // END CELL LOOP ###############################

  // debug
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("air_temp"); vecs.push_back(S_next_->GetFieldData("air_temperature").ptr());
    vnames.push_back("rel hum"); vecs.push_back(S_next_->GetFieldData("relative_humidity").ptr());
    vnames.push_back("Qsw_in"); vecs.push_back(S_next_->GetFieldData("incoming_shortwave_radiation").ptr());
    vnames.push_back("precip_rain"); vecs.push_back(S_next_->GetFieldData("precipitation_rain").ptr());
    vnames.push_back("precip_snow"); vecs.push_back(S_next_->GetFieldData("precipitation_snow").ptr());
    vnames.push_back("T_ground"); vecs.push_back(S_next_->GetFieldData("surface-temperature").ptr());
    vnames.push_back("p_ground"); vecs.push_back(S_next_->GetFieldData("surface-pressure").ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();
    vnames.push_back("snow_ht (old)"); vecs.push_back(S_next_->GetFieldData("snow_depth").ptr());
    vnames.push_back("snow_temp"); vecs.push_back(S_next_->GetFieldData("snow_temperature").ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();
    vnames.push_back("energy_source"); vecs.push_back(S_next_->GetFieldData("surface-conducted_energy_source").ptr());
    vnames.push_back("water_source"); vecs.push_back(S_next_->GetFieldData("surface-mass_source").ptr());
    vnames.push_back("evap flux"); vecs.push_back(S_next_->GetFieldData("evaporative_flux").ptr());
    vnames.push_back("T_water_source"); vecs.push_back(S_next_->GetFieldData("surface-mass_source_temperature").ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();
    vnames.push_back("res (snow_diff)"); vecs.push_back(g->Data().ptr());
    db_->WriteVectors(vnames, vecs, true);
  }

  // mark our other primary variables as changed
  pvfe_esource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wsource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_w_v_source_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wtemp_->SetFieldAsChanged(S_next_.ptr());
  
}


// applies preconditioner to u and returns the result in Pu
int SurfaceBalanceImplicit::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<TreeVector> Pu) {
  *Pu = *u;
  
  return 0;
}


// updates the preconditioner
void
SurfaceBalanceImplicit::UpdatePreconditioner(double t,
        Teuchos::RCP<const TreeVector> up, double h) {}


// error monitor, inf norm is good, this is relative to 1m snow pack
double
SurfaceBalanceImplicit::ErrorNorm(Teuchos::RCP<const TreeVector> u,
        Teuchos::RCP<const TreeVector> du) {
  Teuchos::OSTab tab = vo_->getOSTab();
  double err;
  du->NormInf(&err);
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    *vo_->os() << "ENorm (cells) = " << err << std::endl;
  }
  return err;
}


bool
SurfaceBalanceImplicit::ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
        Teuchos::RCP<TreeVector> u) {
  Epetra_MultiVector& u_vec = *u->Data()->ViewComponent("cell",false);
  unsigned int ncells = u_vec.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    u_vec[0][c] = std::max(0., u_vec[0][c]);
  }

  if (modify_predictor_advance_) {
    Teuchos::RCP<TreeVector> res = Teuchos::rcp(new TreeVector(*u));
    Teuchos::RCP<TreeVector> u0_nc = Teuchos::rcp_const_cast<TreeVector>(u0);
    Functional(S_next_->time()-h, S_next_->time(), u0_nc, u, res);
  }

  return true;
}




} // namespace
} // namespace
