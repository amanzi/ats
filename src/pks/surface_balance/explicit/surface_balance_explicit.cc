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
#include "surface_balance_explicit.hh"

namespace Amanzi {
namespace SurfaceBalance {

SurfaceBalanceExplicit::SurfaceBalanceExplicit(
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           Teuchos::ParameterList& FElist,
           const Teuchos::RCP<TreeVector>& solution) :
    PKPhysicalBase(plist, FElist, solution),
    PKDefaultBase(plist, FElist, solution) {
  // set up additional primary variables -- this is very hacky...
  // -- surface energy source
  Teuchos::ParameterList& esource_sublist =
      FElist.sublist("surface_conducted_energy_source");
  esource_sublist.set("evaluator name", "surface_conducted_energy_source");
  esource_sublist.set("field evaluator type", "primary variable");

  // -- surface mass source
  Teuchos::ParameterList& wsource_sublist =
      FElist.sublist("surface_mass_source");
  wsource_sublist.set("evaluator name", "surface_mass_source");
  wsource_sublist.set("field evaluator type", "primary variable");

  // -- subsurface mass source for VaporFlux at cell center
  Teuchos::ParameterList& w_v_source_sublist =
    FElist.sublist("mass_source");
  w_v_source_sublist.set("evaluator name", "mass_source");
  w_v_source_sublist.set("field evaluator type", "primary variable");

  // -- surface energy temperature
  Teuchos::ParameterList& wtemp_sublist =
      FElist.sublist("surface_mass_source_temperature");
  wtemp_sublist.set("evaluator name", "surface_mass_source_temperature");
  wtemp_sublist.set("field evaluator type", "primary variable");

  // timestep size
  dt_ = plist_->get<double>("max time step", 1.e99);

  // min wind speed
  min_wind_speed_ = plist_->get<double>("minimum wind speed", 1.0);

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
}

// main methods
// -- Setup data.
void
SurfaceBalanceExplicit::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::setup(S);
  subsurf_mesh_ = S->GetMesh(); // needed for VPL, which is treated as subsurface source

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)->
      SetComponent("cell", AmanziMesh::CELL, 1);


  // requirements: other primary variables
  Teuchos::RCP<FieldEvaluator> fm;
  S->RequireField("surface_conducted_energy_source", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_conducted_energy_source");
  fm = S->GetFieldEvaluator("surface_conducted_energy_source");
  pvfe_esource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_esource_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField("surface_mass_source", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_mass_source");
  fm = S->GetFieldEvaluator("surface_mass_source");
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

  S->RequireField("surface_mass_source_temperature", name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("surface_mass_source_temperature");
  fm = S->GetFieldEvaluator("surface_mass_source_temperature");
  pvfe_wtemp_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_wtemp_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceSEB: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  // requirements: independent variables (data from MET)
  S->RequireFieldEvaluator("incoming_shortwave_radiation");
  S->RequireField("incoming_shortwave_radiation")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

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

  // information from ATS data
  S->RequireFieldEvaluator("surface_temperature");
  S->RequireField("surface_temperature")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator("surface_pressure");
  S->RequireField("surface_pressure")->SetMesh(mesh_)
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

  S->RequireFieldEvaluator("surface_porosity");
  S->RequireField("surface_porosity")->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

}

// -- Initialize owned (dependent) variables.
void
SurfaceBalanceExplicit::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBase::initialize(S);

  // initialize snow density
  S->GetFieldData("snow_density",name_)->PutScalar(100.);
  S->GetField("snow_density", name_)->set_initialized();

  // initialize days of no snow
  S->GetFieldData("snow_age",name_)->PutScalar(0.);
  S->GetField("snow_age", name_)->set_initialized();

  // initialize snow temp
  S->GetFieldData("snow_temperature",name_)->PutScalar(0.);
  S->GetField("snow_temperature", name_)->set_initialized();

  // initialize sources, temps
  S->GetFieldData("surface_conducted_energy_source",name_)->PutScalar(0.);
  S->GetField("surface_conducted_energy_source",name_)->set_initialized();

  S->GetFieldData("surface_mass_source",name_)->PutScalar(0.);
  S->GetField("surface_mass_source",name_)->set_initialized();

  S->GetFieldData("mass_source",name_)->PutScalar(0.);
  S->GetField("mass_source",name_)->set_initialized();

  S->GetFieldData("surface_mass_source_temperature",name_)->PutScalar(273.15);
  S->GetField("surface_mass_source_temperature",name_)->set_initialized();
}


// computes the non-linear functional g = g(t,u,udot)
bool
SurfaceBalanceExplicit::advance(double dt) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_inter_->time()
               << " t1 = " << S_next_->time() << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  // pull old snow data
  const Epetra_MultiVector& snow_depth_old = *S_inter_->GetFieldData("snow_depth")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_age_old = *S_inter_->GetFieldData("snow_age")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& snow_dens_old = *S_inter_->GetFieldData("snow_density")
      ->ViewComponent("cell",false);

  // pull current snow data
  Epetra_MultiVector& snow_depth_new = *S_next_->GetFieldData("snow_depth", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& snow_temp_new = *S_next_->GetFieldData("snow_temperature", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& snow_age_new = *S_next_->GetFieldData("snow_age", name_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& snow_dens_new = *S_next_->GetFieldData("snow_density", name_)
      ->ViewComponent("cell",false);

  // pull ATS data
  S_next_->GetFieldEvaluator("surface_temperature")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_temp =
      *S_next_->GetFieldData("surface_temperature")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("surface_pressure")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_pres =
      *S_next_->GetFieldData("surface_pressure")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("ponded_depth")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& ponded_depth =
      *S_next_->GetFieldData("ponded_depth")->ViewComponent("cell", false);

 //  S_next_->GetFieldEvaluator("saturation_liquid")->HasFieldChanged(S_next_.ptr(), name_);
   const Epetra_MultiVector& saturation_liquid =
     *S_next_->GetFieldData("saturation_liquid")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("unfrozen_fraction")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& unfrozen_fraction =
      *S_next_->GetFieldData("unfrozen_fraction")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("surface_porosity")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& surf_porosity =
      *S_next_->GetFieldData("surface_porosity")->ViewComponent("cell", false);


  // pull Met data
  S_next_->GetFieldEvaluator("air_temperature")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& air_temp =
      *S_next_->GetFieldData("air_temperature")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("incoming_shortwave_radiation")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& incoming_shortwave =
      *S_next_->GetFieldData("incoming_shortwave_radiation")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("relative_humidity")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& relative_humidity =
      *S_next_->GetFieldData("relative_humidity")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("wind_speed")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& wind_speed =
      *S_next_->GetFieldData("wind_speed")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("precipitation_rain")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& precip_rain =
      *S_next_->GetFieldData("precipitation_rain")->ViewComponent("cell", false);

  S_next_->GetFieldEvaluator("precipitation_snow")->HasFieldChanged(S_next_.ptr(), name_);
  const Epetra_MultiVector& precip_snow =
      *S_next_->GetFieldData("precipitation_snow")->ViewComponent("cell", false);

  // pull additional primary variable data
  Epetra_MultiVector& surf_energy_flux =
      *S_next_->GetFieldData("surface_conducted_energy_source", name_)->ViewComponent("cell", false);

  Epetra_MultiVector& surf_water_flux =
      *S_next_->GetFieldData("surface_mass_source", name_)->ViewComponent("cell", false);

  Epetra_MultiVector& vapor_flux =
      *S_next_->GetFieldData("mass_source", name_)->ViewComponent("cell", false);
  vapor_flux.PutScalar(0.);

  Epetra_MultiVector& surf_water_flux_temp =
      *S_next_->GetFieldData("surface_mass_source_temperature", name_)->ViewComponent("cell", false);


  unsigned int ncells = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  for (unsigned int c=0; c!=ncells; ++c) {
    if (snow_depth_old[0][c] >= snow_ground_trans_ ||
        snow_depth_old[0][c] < min_snow_trans_) {
      // Evaluate the model as usual
      // Initialize the SEB object
      SEBPhysics::SEB seb;
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
      seb.in.snow_old.ht = snow_depth_old[0][c] < min_snow_trans_ ? 0. : snow_depth_old[0][c];
      seb.in.snow_old.density = snow_dens_old[0][c];
      seb.in.snow_old.age = snow_age_old[0][c];

      seb.out.snow_new.ht = snow_depth_new[0][c];
      seb.out.snow_new.density = snow_dens_new[0][c];
      seb.out.snow_new.age = snow_age_new[0][c];

      seb.in.vp_snow.temp = 273.15;

      // -- snow thermal conductivity parameters
      seb.params.snoK = snow_lambda_fresh_;
      seb.params.snoKX = snow_lambda_exponent_;

      // -- met data
      seb.in.met.Us = std::max(wind_speed[0][c], min_wind_speed_);
      seb.in.met.QswIn = incoming_shortwave[0][c];
      seb.in.met.Ps = precip_snow[0][c];
      seb.in.met.Pr = precip_rain[0][c];
      seb.in.met.vp_air.temp = air_temp[0][c];
      seb.in.met.vp_air.relative_humidity = relative_humidity[0][c];

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
      SEBPhysics::CalculateSurfaceBalance(seb, true);

      // Pull the output
      // -- fluxes
      surf_energy_flux[0][c] = seb.out.eb.fQc;
      surf_water_flux[0][c] = seb.out.mb.MWg;
      surf_water_flux_temp[0][c] = seb.out.mb.MWg_temp;

      // -- vapor flux to cells
      //     surface vapor flux is treated as a volumetric source for the subsurface.
//      AmanziMesh::Entity_ID subsurf_f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
//      AmanziMesh::Entity_ID_List cells;
      subsurf_mesh_->face_get_cells(subsurf_f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);
      // surface mass sources are in m^3 water / (m^2 s)
      // subsurface mass sources are in mol water / (m^3 s)
      vapor_flux[0][cells[0]] = seb.out.mb.MWg_subsurf
          * mesh_->cell_volume(c) * seb.in.vp_ground.density_w / 0.0180153
          / subsurf_mesh_->cell_volume(cells[0]);

      // -- snow properties
      snow_depth_new[0][c] = seb.out.snow_new.ht;
      snow_age_new[0][c] = seb.out.snow_new.age;
      snow_dens_new[0][c] = seb.out.snow_new.density;
      snow_temp_new[0][c] = seb.in.vp_snow.temp;


    } else {
      // Evaluate the model twice -- once as bare ground, once with snow, using
      // an area-averaged subgrid model to smooth between the two end-members.
      // The area-weighting parameter is theta:
      double theta = snow_depth_old[0][c] / snow_ground_trans_;

      // Evaluate the model as usual
      // Initialize the SEB object
      SEBPhysics::SEB seb;
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
      seb.in.snow_old.density = snow_dens_old[0][c];
      seb.in.snow_old.age = snow_age_old[0][c];

      seb.out.snow_new.ht = snow_ground_trans_;
      seb.out.snow_new.density = snow_dens_new[0][c];
      seb.out.snow_new.age = snow_age_new[0][c];

      seb.in.vp_snow.temp = 273.15;

      // -- met data
      seb.in.met.Us = std::max(wind_speed[0][c], min_wind_speed_);
      seb.in.met.QswIn = incoming_shortwave[0][c];
      seb.in.met.Ps = precip_snow[0][c];
      seb.in.met.Pr = precip_rain[0][c];
      seb.in.met.vp_air.temp = air_temp[0][c];
      seb.in.met.vp_air.relative_humidity = relative_humidity[0][c];

      // -- snow thermal conductivity parameters
      seb.params.snoK = snow_lambda_fresh_;
      seb.params.snoKX = snow_lambda_exponent_;

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
      SEBPhysics::CalculateSurfaceBalance(seb, true);
      SEBPhysics::CalculateSurfaceBalance(seb_bare, true);

      // Pull the output
      // -- fluxes
      surf_energy_flux[0][c] = theta * seb.out.eb.fQc + (1-theta) * seb_bare.out.eb.fQc;
      surf_water_flux[0][c] = theta * seb.out.mb.MWg + (1-theta) * seb_bare.out.mb.MWg;
      if (std::abs(surf_water_flux[0][c]) > 0.) {
        surf_water_flux_temp[0][c] = (theta * seb.out.mb.MWg * seb.out.mb.MWg_temp
                                      + (1-theta) * seb_bare.out.mb.MWg * seb_bare.out.mb.MWg_temp) / surf_water_flux[0][c];
      } else {
        surf_water_flux_temp[0][c] = seb.out.mb.MWg_temp;
      }

      // -- vapor flux to cells
      //     surface vapor flux is treated as a volumetric source for the subsurface.
//      AmanziMesh::Entity_ID subsurf_f = mesh_->entity_get_parent(AmanziMesh::CELL, c);
//      AmanziMesh::Entity_ID_List cells;
      subsurf_mesh_->face_get_cells(subsurf_f, AmanziMesh::OWNED, &cells);
      ASSERT(cells.size() == 1);
      // surface mass sources are in m^3 water / (m^2 s)
      // subsurface mass sources are in mol water / (m^3 s)
      double mean_flux = theta * seb.out.mb.MWg_subsurf + (1-theta) * seb_bare.out.mb.MWg_subsurf;
      vapor_flux[0][cells[0]] = mean_flux
          * mesh_->cell_volume(c) * seb.in.vp_ground.density_w / 0.0180153
          / subsurf_mesh_->cell_volume(cells[0]);

      // -- snow properties: SWE averaged
      double total_swe = theta * seb.out.snow_new.ht * seb.out.snow_new.density / seb.in.vp_ground.density_w
          + (1-theta) * seb_bare.out.snow_new.ht * seb_bare.out.snow_new.density / seb_bare.in.vp_ground.density_w;
      snow_depth_new[0][c] = theta * seb.out.snow_new.ht + (1-theta) * seb_bare.out.snow_new.ht;
      if (snow_depth_new[0][c] > 0.) {
        snow_age_new[0][c] = (theta * seb.out.snow_new.ht * seb.out.snow_new.age
                              + (1-theta) * seb_bare.out.snow_new.ht * seb_bare.out.snow_new.age) / snow_depth_new[0][c];
        snow_dens_new[0][c] = total_swe * seb.in.vp_ground.density_w / snow_depth_new[0][c];
      } else {
        snow_age_new[0][c] = 0.;
        snow_dens_new[0][0] = seb.out.snow_new.density;
      }
      snow_temp_new[0][c] = seb.in.vp_snow.temp;
    }
  }

  // Mark primary variables as changed.
  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
  pvfe_esource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wsource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_w_v_source_->SetFieldAsChanged(S_next_.ptr());
  pvfe_wtemp_->SetFieldAsChanged(S_next_.ptr());

  // debug
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
    vnames.push_back("air_temp"); vecs.push_back(S_next_->GetFieldData("air_temperature").ptr());
    vnames.push_back("Qsw_in"); vecs.push_back(S_next_->GetFieldData("incoming_shortwave_radiation").ptr());
    vnames.push_back("precip_rain"); vecs.push_back(S_next_->GetFieldData("precipitation_rain").ptr());
    vnames.push_back("precip_snow"); vecs.push_back(S_next_->GetFieldData("precipitation_snow").ptr());
    vnames.push_back("T_ground"); vecs.push_back(S_next_->GetFieldData("surface_temperature").ptr());
    vnames.push_back("p_ground"); vecs.push_back(S_next_->GetFieldData("surface_pressure").ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();
    vnames.push_back("snow_ht"); vecs.push_back(S_next_->GetFieldData("snow_depth").ptr());
    vnames.push_back("snow_temp"); vecs.push_back(S_next_->GetFieldData("snow_temperature").ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();
    vnames.push_back("energy_source"); vecs.push_back(S_next_->GetFieldData("surface_conducted_energy_source").ptr());
    vnames.push_back("water_source"); vecs.push_back(S_next_->GetFieldData("surface_mass_source").ptr());
    vnames.push_back("surface_vapor_source"); vecs.push_back(S_next_->GetFieldData("mass_source").ptr());
    vnames.push_back("T_water_source"); vecs.push_back(S_next_->GetFieldData("surface_mass_source_temperature").ptr());
    db_->WriteVectors(vnames, vecs, true);
  }


  return false;
}

} // namespace
} // namespace
