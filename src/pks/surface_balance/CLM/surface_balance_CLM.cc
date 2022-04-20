/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Uses CLM for determining transpiration, evaporation, snowmelt, etc.

/*!

.. note: This is currently a PK.  But it acts like an evaluator.  Much of this
    code should get moved into an evaluator.

Based on the Colorado/Common/Community Land Model, an old variant of CLM that
has been updated and maintained by the ParFlow group.

CLM provides all surface processes, including snowpack evolution, energy and
water sources, etc.

.. note: Currently this is a water-only model -- it internally does its own
    energy equation.  One could refactor CLM to split out this energy balance
    as well, allowing us to use this with ATS's energy equations, but that is
    currently not possible.


 */

#include <cmath>

#include "ats_clm_interface.hh"
#include "surface_balance_CLM.hh"


namespace Amanzi {
namespace SurfaceBalance {

SurfaceBalanceCLM::SurfaceBalanceCLM(Teuchos::ParameterList& pk_tree,
                                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                     const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& solution):
  PK(pk_tree, global_list,  S, solution),
  PK_Physical_Default(pk_tree, global_list,  S, solution),
  my_next_time_(0.)
{
  domain_ss_ = Keys::readDomainHint(*plist_, domain_, "surface", "subsurface");

  // set up primary variables for surface/subsurface sources.  CLM keeps its
  // own internal state, violating most ATS principles, but for now we'll hack it in.
  // -- surface water sources
  Teuchos::ParameterList& wsource_sublist =
    S->GetEvaluatorList(Keys::getKey(domain_,"water_source"));
  wsource_sublist.set("evaluator name", Keys::getKey(domain_,"water_source"));
  wsource_sublist.set("field evaluator type", "primary variable");

  // -- subsurface water source transpiration
  Teuchos::ParameterList& w_v_source_sublist =
    S->GetEvaluatorList(Keys::getKey(domain_ss_,"water_source"));
  w_v_source_sublist.set("evaluator name", Keys::getKey(domain_ss_,"water_source"));
  w_v_source_sublist.set("field evaluator type", "primary variable");

  // CLM timestep
  dt_ = plist_->get<double>("time step size [s]");
}

// main methods
// -- Setup data.
void
SurfaceBalanceCLM::Setup(const Teuchos::Ptr<State>& S) {
  PK_Physical_Default::Setup(S);
  subsurf_mesh_ = S->GetMesh(domain_ss_); // needed for VPL, which is treated as subsurface source

  // requirements: primary variable
  S->RequireField(key_, name_)->SetMesh(mesh_)->
      SetComponent("cell", AmanziMesh::CELL, 1);

  // requirements: other primary variables
  Teuchos::RCP<FieldEvaluator> fm;
  S->RequireField(Keys::getKey(domain_,"water_source"), name_)->SetMesh(mesh_)
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_,"water_source"));
  fm = S->GetFieldEvaluator(Keys::getKey(domain_,"water_source"));
  pvfe_wsource_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_wsource_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceCLM: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }

  S->RequireField(Keys::getKey(domain_ss_,"water_source"), name_)->SetMesh(subsurf_mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_ss_,"water_source"));
  fm = S->GetFieldEvaluator(Keys::getKey(domain_ss_,"water_source"));
  pvfe_w_sub_source_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  if (pvfe_w_sub_source_ == Teuchos::null) {
    Errors::Message message("SurfaceBalanceCLM: error, failure to initialize primary variable");
    Exceptions::amanzi_throw(message);
  }
  
  // requirements: energy balance diagnostic variables
  S->RequireField(Keys::getKey(domain_,"evaporative_flux"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"evaporative_flux"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"qE_latent_heat"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"qE_latent_heat"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"qE_sensible_heat"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"qE_sensible_heat"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"qE_lw_out"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"qE_lw_out"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"qE_conducted_soil"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"qE_conducted_soil"),name_)->set_io_checkpoint(false);

  // requirements: other diagnostics
  S->RequireField(Keys::getKey(domain_,"snow_swe"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"snow_swe"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"canopy_storage"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"canopy_storage"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"temperature_skin"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"temperature_skin"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_,"temperature_leaf"),name_)->SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_,"temperature_leaf"),name_)->set_io_checkpoint(false);
  S->RequireField(Keys::getKey(domain_ss_,"temperature_soil"),name_)->SetMesh(subsurf_mesh_)
    ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->GetField(Keys::getKey(domain_ss_,"temperature_soil"),name_)->set_io_checkpoint(false);

  
  // requirements: independent variables (data from MET)
  S->RequireFieldEvaluator(Keys::getKey(domain_, "incoming_shortwave_radiation"));
  S->RequireField(Keys::getKey(domain_, "incoming_shortwave_radiation"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"incoming_longwave_radiation"));
  S->RequireField(Keys::getKey(domain_,"incoming_longwave_radiation"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"air_temperature"));
  S->RequireField(Keys::getKey(domain_,"air_temperature"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"relative_humidity"));
  S->RequireField(Keys::getKey(domain_,"relative_humidity"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"wind_speed"));
  S->RequireField(Keys::getKey(domain_,"wind_speed"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"precipitation_rain"));
  S->RequireField(Keys::getKey(domain_,"precipitation_rain"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_,"precipitation_snow"));
  S->RequireField(Keys::getKey(domain_,"precipitation_snow"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirements: soil state
  S->RequireFieldEvaluator(Keys::getKey(domain_,"pressure"));
  S->RequireField(Keys::getKey(domain_,"pressure"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_ss_,"porosity"));
  S->RequireField(Keys::getKey(domain_ss_,"porosity"))->SetMesh(subsurf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  S->RequireFieldEvaluator(Keys::getKey(domain_ss_,"saturation_liquid"));
  S->RequireField(Keys::getKey(domain_ss_,"saturation_liquid"))->SetMesh(subsurf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

  // requirements: soil properties
  S->RequireFieldEvaluator(Keys::getKey(domain_ss_, "sand_fraction"));
  S->RequireField(Keys::getKey(domain_ss_,"sand_fraction"))->SetMesh(subsurf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_ss_, "clay_fraction"));
  S->RequireField(Keys::getKey(domain_ss_,"clay_fraction"))->SetMesh(subsurf_mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_, "color_index"));
  S->RequireField(Keys::getKey(domain_,"color_index"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(Keys::getKey(domain_, "pft_index"));
  S->RequireField(Keys::getKey(domain_,"pft_index"))->SetMesh(mesh_)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  
  
  // Set up the CLM object
  ATS::CLM::init(subsurf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED),
                 mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED),
                 2, mesh_->get_comm()->MyPID(), 3);
}

// -- Initialize owned (dependent) variables.
void
SurfaceBalanceCLM::Initialize(const Teuchos::Ptr<State>& S) {
  PK_Physical_Default::Initialize(S);

  Teuchos::ParameterList& ic_list = plist_->sublist("initial condition");
  double snow_depth = ic_list.get<double>("initial snow depth [m]");
  double temp = ic_list.get<double>("initial temperature [K]");
  double year = ic_list.get<double>("initial time [yr]");
  ATS::CLM::set_zero_time(year);
  ATS::CLM::set_initial_state(temp, snow_depth);

  // lat/lon
  auto latlon = plist_->get<Teuchos::Array<double>>("latitude,longitude [degrees]");
  int ncols = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  double latlon_arr[ncols][2];
  for (int i=0; i!=ncols; ++i) {
    latlon_arr[i][0] = latlon[0];
    latlon_arr[i][1] = latlon[1];
  }
  S->GetFieldEvaluator(Keys::getKey(domain_ss_, "sand_fraction"))
      ->HasFieldChanged(S.ptr(), name_);
  auto& sand = *S->GetFieldData(Keys::getKey(domain_ss_, "sand_fraction"))
               ->ViewComponent("cell", false);
  S->GetFieldEvaluator(Keys::getKey(domain_ss_, "clay_fraction"))
      ->HasFieldChanged(S.ptr(), name_);
  auto& clay = *S->GetFieldData(Keys::getKey(domain_ss_, "clay_fraction"))
               ->ViewComponent("cell", false);
  S->GetFieldEvaluator(Keys::getKey(domain_, "color_index"))
      ->HasFieldChanged(S.ptr(), name_);
  auto& color = *S->GetFieldData(Keys::getKey(domain_, "color_index"))
               ->ViewComponent("cell", false);
  std::vector<int> color_index(ncols);
  AMANZI_ASSERT(color.MyLength() == ncols);
  for (int i=0; i!=ncols; ++i) color_index[i] = std::round(color[0][i]);

  // pft tile
  S->GetFieldEvaluator(Keys::getKey(domain_, "pft_index"))
      ->HasFieldChanged(S.ptr(), name_);
  auto& pft = *S->GetFieldData(Keys::getKey(domain_, "pft_index"))
               ->ViewComponent("cell", false);
  AMANZI_ASSERT(pft.MyLength() == ncols);
  double fractional_ground[ncols][NUM_LC_CLASSES];
  for (int i=0; i!=ncols; ++i) {
    for (int j=0; j!=NUM_LC_CLASSES; ++j) {
      fractional_ground[i][j] = j == std::round(pft[0][i]) ? 1. : 0.;
    }
  }
  ATS::CLM::set_ground_properties(&latlon_arr[0][0], sand, clay, color_index, &fractional_ground[0][0]);

  ATS::CLM::setup_begin();
  Epetra_MultiVector dz(subsurf_mesh_->cell_map(false), 1);
  for (int col=0; col!=ncols; ++col) {
    auto& faces = subsurf_mesh_->faces_of_column(col);
    auto& cells = subsurf_mesh_->cells_of_column(col);
    for (int i=0; i!=cells.size(); ++i) {
      dz[0][cells[i]] = subsurf_mesh_->face_centroid(faces[i])[2] -
                     subsurf_mesh_->face_centroid(faces[i+1])[2];
      AMANZI_ASSERT(dz[0][cells[i]] > 0.);
    }    
  }
  ATS::CLM::set_dz(dz);
  ATS::CLM::set_et_controls(1, 2, 0.1, 1.0, 0.1);
  ATS::CLM::setup_end();
  ATS::CLM::set_dz(dz);

  // set as intialized the sources
  S->GetFieldData(Keys::getKey(domain_,"water_source"),name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_,"water_source"),name_)->set_initialized();
  S->GetFieldData(Keys::getKey(domain_ss_,"water_source"),name_)->PutScalar(0.);
  S->GetField(Keys::getKey(domain_ss_,"water_source"),name_)->set_initialized();
  S->GetFieldData(Keys::getKey(domain_,"snow_depth"),name_)->PutScalar(snow_depth);
  S->GetField(Keys::getKey(domain_,"snow_depth"),name_)->set_initialized();

  // set as intialized the diagnostics
  S->GetField(Keys::getKey(domain_,"evaporative_flux"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"qE_latent_heat"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"qE_sensible_heat"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"qE_lw_out"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"qE_conducted_soil"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"snow_swe"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"canopy_storage"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"temperature_skin"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_,"temperature_leaf"),name_)->set_initialized();
  S->GetField(Keys::getKey(domain_ss_,"temperature_soil"),name_)->set_initialized();
}


bool
SurfaceBalanceCLM::AdvanceStep(double t_old, double t_new, bool reinit) {
  if (t_new <= my_next_time_) {
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
                 << "BIG STEP still good!" << std::endl
                 << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;
    return false;
  }
  
  Teuchos::OSTab tab = vo_->getOSTab();

  bool debug = false;
  Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
  int rank = mesh_->get_comm()->MyPID();

  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl
               << "BIG STEP Advancing: t0 = " << t_old
               << " t1 = " << t_old + dt_ << " h = " << dt_ << std::endl
               << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << std::endl;


  // Set the state
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_ss_, "pressure"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& pressure = *S_inter_->GetFieldData(Keys::getKey(domain_ss_, "pressure"))
                                       ->ViewComponent("cell", false);
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_ss_, "pressure"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& porosity = *S_inter_->GetFieldData(Keys::getKey(domain_ss_, "porosity"))
                                       ->ViewComponent("cell", false);
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_ss_, "saturation_liquid"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& sl = *S_inter_->GetFieldData(Keys::getKey(domain_ss_, "saturation_liquid"))
                                       ->ViewComponent("cell", false);
  double patm = *S_inter_->GetScalarData("atmospheric_pressure");

  ATS::CLM::set_wc(porosity, sl);
  ATS::CLM::set_tksat_from_porosity(porosity);
  ATS::CLM::set_pressure(pressure, patm);

  // set the forcing
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_, "incoming_shortwave_radiation"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& qSW = *S_inter_->GetFieldData(Keys::getKey(domain_,
          "incoming_shortwave_radiation"))->ViewComponent("cell", false);
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_, "incoming_longwave_radiation"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& qLW = *S_inter_->GetFieldData(Keys::getKey(domain_,
          "incoming_longwave_radiation"))->ViewComponent("cell", false);
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_, "precipitation_snow"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& pSnow = *S_inter_->GetFieldData(Keys::getKey(domain_,
          "precipitation_snow"))->ViewComponent("cell", false);
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_, "precipitation_rain"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& pRain = *S_inter_->GetFieldData(Keys::getKey(domain_,
          "precipitation_rain"))->ViewComponent("cell", false);
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_, "air_temperature"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& air_temp = *S_inter_->GetFieldData(Keys::getKey(domain_,
          "air_temperature"))->ViewComponent("cell", false);
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_, "relative_humidity"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& rel_hum = *S_inter_->GetFieldData(Keys::getKey(domain_,
          "relative_humidity"))->ViewComponent("cell", false);
  S_inter_->GetFieldEvaluator(Keys::getKey(domain_, "wind_speed"))
      ->HasFieldChanged(S_inter_.ptr(), name_);
  const Epetra_MultiVector& wind_speed = *S_inter_->GetFieldData(Keys::getKey(domain_,
          "wind_speed"))->ViewComponent("cell", false);
  ATS::CLM::set_met_data(qSW, qLW, pRain, pSnow, air_temp, rel_hum, wind_speed, patm);

  // set the start time, endtime
  ATS::CLM::advance_time(S_inter_->cycle(), t_old, dt_); //units in seconds

  // get diagnostics
  Epetra_MultiVector& latent_heat = *S_next_->GetFieldData(Keys::getKey(domain_,
          "qE_latent_heat"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& sensible_heat = *S_next_->GetFieldData(Keys::getKey(domain_,
          "qE_sensible_heat"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& lw_out = *S_next_->GetFieldData(Keys::getKey(domain_,
          "qE_lw_out"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& cond = *S_next_->GetFieldData(Keys::getKey(domain_,
          "qE_conducted_soil"), name_)->ViewComponent("cell", false);
  ATS::CLM::get_ground_energy_fluxes(latent_heat,sensible_heat,lw_out, cond);

  // more diagnostics
  Epetra_MultiVector& swe = *S_next_->GetFieldData(Keys::getKey(domain_,
          "snow_swe"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& snow_depth = *S_next_->GetFieldData(key_, name_)
                                   ->ViewComponent("cell", false);
  Epetra_MultiVector& canopy_storage = *S_next_->GetFieldData(Keys::getKey(domain_,
          "canopy_storage"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& temperature_skin = *S_next_->GetFieldData(Keys::getKey(domain_,
          "temperature_skin"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& temperature_leaf = *S_next_->GetFieldData(Keys::getKey(domain_,
          "temperature_leaf"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& temperature_soil = *S_next_->GetFieldData(Keys::getKey(domain_ss_,
          "temperature_soil"), name_)->ViewComponent("cell", false);
  ATS::CLM::get_diagnostics(swe, snow_depth, canopy_storage, temperature_skin,
                            temperature_leaf, temperature_soil);
  
  
  // get output
  Epetra_MultiVector& surf_source = *S_next_->GetFieldData(Keys::getKey(domain_,
          "water_source"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& sub_source = *S_next_->GetFieldData(Keys::getKey(domain_ss_,
          "water_source"), name_)->ViewComponent("cell", false);
  ATS::CLM::get_total_mass_fluxes(surf_source, sub_source);

  Epetra_MultiVector& surf_source_old = *S_inter_->GetFieldData(Keys::getKey(domain_,
          "water_source"), name_)->ViewComponent("cell", false);
  Epetra_MultiVector& sub_source_old = *S_inter_->GetFieldData(Keys::getKey(domain_ss_,
          "water_source"), name_)->ViewComponent("cell", false);
  surf_source_old = surf_source;
  sub_source_old = sub_source;

  pvfe_wsource_->SetFieldAsChanged(S_next_.ptr());
  pvfe_w_sub_source_->SetFieldAsChanged(S_next_.ptr());
  my_next_time_ = t_old + dt_;

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector< Teuchos::Ptr<const CompositeVector> > vecs;

    vnames.push_back("inc shortwave radiation [W/m^2]");
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_, "incoming_shortwave_radiation")).ptr());
    vnames.push_back("inc longwave radiation [W/m^2]");
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_, "incoming_longwave_radiation")).ptr());
    vnames.push_back("inc latent heat [W/m^2]"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"qE_latent_heat")).ptr());
    vnames.push_back("inc sensible heat [W/m^2]"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"qE_sensible_heat")).ptr());
    vnames.push_back("out longwave radiation [W/m^2]"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"qE_lw_out")).ptr());
    vnames.push_back("out conducted soil [W/m^2]"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"qE_conducted_soil")).ptr());

    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();
    
    vnames.push_back("surface water source [m/s]"); 
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"water_source")).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("snow depth [m]");
    vecs.push_back(S_next_->GetFieldData(key_).ptr());
    vnames.push_back("snow swe [m]");
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_, "snow_swe")).ptr());
    vnames.push_back("canopy storage [m]");
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_, "canopy_storage")).ptr());
    vnames.push_back("skin temperature [K]");
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_, "temperature_skin")).ptr());
    vnames.push_back("leaf temperature [K]");
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_, "temperature_leaf")).ptr());
    
  }

  return false;
}


} // namespace
} // namespace
