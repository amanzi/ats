/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Uses CLM for determining transpiration, evaporation, snowmelt, etc.
/*!

.. note: This is currently a PK.  But it acts like an evaluator.  Much of this
    code should get moved into an evaluator.

Based on the Community Land Model, an old variant of CLM that
has been updated and maintained by the ParFlow group.

CLM provides all surface processes, including snowpack evolution, energy and
water sources, etc.

.. note: Currently this is a water-only model -- it internally does its own
    energy equation.  One could refactor CLM to split out this energy balance
    as well, allowing us to use this with ATS's energy equations, but that is
    currently not possible.

*/

#include <cmath>

#include "pk_helpers.hh"
#include "ats_clm_interface.hh"
#include "surface_balance_CLM.hh"

namespace Amanzi {
namespace SurfaceBalance {

SurfaceBalanceCLM::SurfaceBalanceCLM(Teuchos::ParameterList& pk_tree,
                                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                     const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, global_list, S, solution), PK_Physical_Default(pk_tree, global_list, S, solution)
{
  domain_ss_ = Keys::readDomainHint(*plist_, domain_, "surface", "subsurface");
  domain_snow_ = Keys::readDomainHint(*plist_, domain_, "surface", "snow");
  domain_can_ = Keys::readDomainHint(*plist_, domain_, "surface", "canopy");

  // primary variables
  key_ = Keys::readKey(*plist_, domain_snow_, "snow depth", key_);
  surf_water_src_key_ = Keys::readKey(*plist_, domain_, "surface water source", "water_source");
  ss_water_src_key_ = Keys::readKey(*plist_, domain_ss_, "subsurface water source", "water_source");

  // diagnostic keys
  // evap_flux_key_ = Keys::readKey(*plist_, domain_, "evaporative flux", "evaporative_flux");
  qE_lh_key_ = Keys::readKey(*plist_, domain_, "latent heat of evaporation", "qE_latent_heat");
  qE_sh_key_ = Keys::readKey(*plist_, domain_, "sensible heat flux", "qE_sensible_heat");
  qE_lw_out_key_ = Keys::readKey(*plist_, domain_, "outgoing longwave radiation", "qE_lw_out");
  qE_cond_key_ = Keys::readKey(*plist_, domain_, "conducted energy flux", "qE_conducted");

  snow_swe_key_ = Keys::readKey(*plist_, domain_snow_, "snow water equivalent", "water_equivalent");
  can_wc_key_ = Keys::readKey(*plist_, domain_can_, "canopy water content", "water_content");
  surf_temp_key_ = Keys::readKey(*plist_, domain_, "surface temperature", "temperature");
  soil_temp_key_ = Keys::readKey(*plist_, domain_ss_, "soil temperature", "temperature");
  can_temp_key_ = Keys::readKey(*plist_, domain_can_, "canopy temperature", "temperature");

  // dependencies
  // -- met data
  met_sw_key_ =
    Keys::readKey(*plist_, domain_, "incoming shortwave radiation", "incoming_shortwave_radiation");
  met_lw_key_ =
    Keys::readKey(*plist_, domain_, "incoming longwave radiation", "incoming_longwave_radiation");
  met_air_temp_key_ = Keys::readKey(*plist_, domain_, "air temperature", "air_temperature");
  met_vp_air_key_ = Keys::readKey(*plist_, domain_, "vapor pressure air", "vapor_pressure_air");
  met_wind_speed_key_ = Keys::readKey(*plist_, domain_, "wind speed", "wind_speed");
  met_prain_key_ = Keys::readKey(*plist_, domain_, "precipitation rain", "precipitation_rain");
  met_psnow_key_ = Keys::readKey(*plist_, domain_snow_, "precipitation snow", "precipitation");

  // soil state
  pres_key_ = Keys::readKey(*plist_, domain_ss_, "pressure", "pressure");
  sl_key_ = Keys::readKey(*plist_, domain_ss_, "saturation_liquid", "saturation_liquid");

  // soil properties
  sand_frac_key_ = Keys::readKey(*plist_, domain_ss_, "sand fraction", "sand_fraction");
  silt_frac_key_ = Keys::readKey(*plist_, domain_ss_, "silt fraction", "silt_fraction");
  clay_frac_key_ = Keys::readKey(*plist_, domain_ss_, "clay fraction", "clay_fraction");
  poro_key_ = Keys::readKey(*plist_, domain_ss_, "porosity", "porosity");

  // surface properties
  color_index_key_ = Keys::readKey(*plist_, domain_ss_, "color index", "color_index");
  pft_index_key_ = Keys::readKey(*plist_, domain_ss_, "PFT index", "pft_index");

  // CLM timestep
  dt_ = plist_->get<double>("time step size [s]");
}

// main methods
// -- Setup data.
void
SurfaceBalanceCLM::Setup()
{
  PK_Physical_Default::Setup();
  auto subsurf_mesh = S_->GetMesh(domain_ss_);

  // requirements: primary variable
  // -- snow depth
  S_->Require<CompositeVector, CompositeVectorSpace>(key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // requirements: other primary variables
  // -- surface water source  -- note we keep old and new in case of Crank-Nicholson Richards PK
  S_->Require<CompositeVector, CompositeVectorSpace>(surf_water_src_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireEvaluatorPrimary(surf_water_src_key_, tag_next_, *S_);

  // -- subsurface water source  -- note we keep old and new in case of Crank-Nicholson Richards PK
  S_->Require<CompositeVector, CompositeVectorSpace>(ss_water_src_key_, tag_next_, name_)
    .SetMesh(subsurf_mesh)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  requireEvaluatorPrimary(ss_water_src_key_, tag_next_, *S_);

  // set requirements on dependencies
  SetupDependencies_(tag_next_);

  // Set up the CLM object
  ATS::CLM::init(subsurf_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED),
                 mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED),
                 2,
                 mesh_->getComm()->MyPID(),
                 3);
}

void
SurfaceBalanceCLM::SetupDependencies_(const Tag& tag)
{
  auto subsurf_mesh = S_->GetMesh(domain_ss_);
  auto snow_mesh = S_->GetMesh(domain_snow_);
  auto can_mesh = S_->GetMesh(domain_can_);

  // requirements: energy balance diagnostic variables.  Only at the new time.
  // No evaluators for now?
  // S_->Require<CompositeVector,CompositeVectorSpace>(evap_flux_key_, tag_next_, name_)
  //   .SetMesh(mesh_)->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  // S_->GetRecord(evap_flux_key_, tag_next_).set_io_checkpoint(false);

  S_->Require<CompositeVector, CompositeVectorSpace>(qE_lh_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->GetRecordW(qE_lh_key_, tag_next_, name_).set_io_checkpoint(false);

  S_->Require<CompositeVector, CompositeVectorSpace>(qE_sh_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->GetRecordW(qE_sh_key_, tag_next_, name_).set_io_checkpoint(false);

  S_->Require<CompositeVector, CompositeVectorSpace>(qE_lw_out_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->GetRecordW(qE_lw_out_key_, tag_next_, name_).set_io_checkpoint(false);

  S_->Require<CompositeVector, CompositeVectorSpace>(qE_cond_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->GetRecordW(qE_cond_key_, tag_next_, name_).set_io_checkpoint(false);

  // requirements: other diagnostics
  S_->Require<CompositeVector, CompositeVectorSpace>(snow_swe_key_, tag_next_, name_)
    .SetMesh(snow_mesh)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->GetRecordW(snow_swe_key_, tag_next_, name_).set_io_checkpoint(false);

  S_->Require<CompositeVector, CompositeVectorSpace>(can_wc_key_, tag_next_, name_)
    .SetMesh(can_mesh)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->GetRecordW(can_wc_key_, tag_next_, name_).set_io_checkpoint(false);

  S_->Require<CompositeVector, CompositeVectorSpace>(surf_temp_key_, tag_next_, name_)
    .SetMesh(mesh_)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->GetRecordW(surf_temp_key_, tag_next_, name_).set_io_checkpoint(false);

  S_->Require<CompositeVector, CompositeVectorSpace>(soil_temp_key_, tag_next_, name_)
    .SetMesh(subsurf_mesh)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->GetRecordW(soil_temp_key_, tag_next_, name_).set_io_checkpoint(false);

  S_->Require<CompositeVector, CompositeVectorSpace>(can_temp_key_, tag_next_, name_)
    .SetMesh(can_mesh)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->GetRecordW(can_temp_key_, tag_next_, name_).set_io_checkpoint(false);

  // requirements: independent variables (data from MET)
  S_->RequireEvaluator(met_sw_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(met_sw_key_, tag)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator(met_lw_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(met_lw_key_, tag)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator(met_air_temp_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(met_air_temp_key_, tag)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator(met_vp_air_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(met_vp_air_key_, tag)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator(met_wind_speed_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(met_wind_speed_key_, tag)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator(met_prain_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(met_prain_key_, tag)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S_->RequireEvaluator(met_psnow_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(met_psnow_key_, tag)
    .SetMesh(snow_mesh)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // requirements: soil state
  S_->RequireEvaluator(pres_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(pres_key_, tag)
    .SetMesh(subsurf_mesh)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(sl_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(sl_key_, tag)
    .SetMesh(subsurf_mesh)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // requirements: soil properties
  S_->RequireEvaluator(poro_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(poro_key_, tag)
    .SetMesh(subsurf_mesh)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(sand_frac_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(sand_frac_key_, tag)
    .SetMesh(subsurf_mesh)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(silt_frac_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(silt_frac_key_, tag)
    .SetMesh(subsurf_mesh)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(clay_frac_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(clay_frac_key_, tag)
    .SetMesh(subsurf_mesh)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // requirements: surface properties
  S_->RequireEvaluator(color_index_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(color_index_key_, tag)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(pft_index_key_, tag);
  S_->Require<CompositeVector, CompositeVectorSpace>(pft_index_key_, tag)
    .SetMesh(mesh_)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}


// -- Initialize owned (dependent) variables.
void
SurfaceBalanceCLM::Initialize()
{
  PK_Physical_Default::Initialize();
  InitializeCLM_(tag_next_);
  InitializePrimaryVariables_(tag_next_);
}


void
SurfaceBalanceCLM::InitializeCLM_(const Tag& tag)
{
  // Initialize the CLM instance
  Teuchos::ParameterList& ic_list = plist_->sublist("initial condition");
  double snow_depth = ic_list.get<double>("initial snow depth [m]");
  double temp = ic_list.get<double>("initial temperature [K]");
  double year = ic_list.get<double>("initial time [yr]");
  ATS::CLM::set_zero_time(year);
  ATS::CLM::set_initial_state(temp, snow_depth);

  // lat/lon
  auto latlon = plist_->get<Teuchos::Array<double>>("latitude,longitude [degrees]");
  int ncols = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
  double latlon_arr[ncols][2];
  for (int i = 0; i != ncols; ++i) {
    latlon_arr[i][0] = latlon[0];
    latlon_arr[i][1] = latlon[1];
  }

  // soil properties
  S_->GetEvaluator(sand_frac_key_, tag).Update(*S_, name_);
  auto& sand = *S_->Get<CompositeVector>(sand_frac_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(silt_frac_key_, tag).Update(*S_, name_);
  auto& silt = *S_->Get<CompositeVector>(silt_frac_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(clay_frac_key_, tag).Update(*S_, name_);
  auto& clay = *S_->Get<CompositeVector>(clay_frac_key_, tag).ViewComponent("cell", false);

  S_->GetEvaluator(color_index_key_, tag).Update(*S_, name_);
  auto& color_index_tmp =
    *S_->Get<CompositeVector>(color_index_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(pft_index_key_, tag).Update(*S_, name_);
  auto& pft_index_tmp = *S_->Get<CompositeVector>(pft_index_key_, tag).ViewComponent("cell", false);

  std::vector<int> color_index(ncols);
  for (int i = 0; i != ncols; ++i) color_index[i] = std::round(color_index_tmp[0][i]);

  double pft_fraction[ncols][NUM_LC_CLASSES];
  for (int i = 0; i != ncols; ++i) {
    for (int j = 0; j != NUM_LC_CLASSES; ++j) {
      pft_fraction[i][j] = j == std::round(pft_index_tmp[0][i]) ? 1. : 0.;
    }
  }
  ATS::CLM::set_ground_properties(&latlon_arr[0][0], sand, clay, color_index, &pft_fraction[0][0]);

  // CLM setup stage
  ATS::CLM::setup_begin();
  auto subsurf_mesh = S_->GetMesh(domain_ss_);
  Epetra_MultiVector dz(subsurf_mesh->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
  for (int col = 0; col != ncols; ++col) {
    auto& faces = subsurf_mesh->columns.getFaces(col);
    auto& cells = subsurf_mesh->columns.getCells(col);
    for (int i = 0; i != cells.size(); ++i) {
      dz[0][cells[i]] =
        subsurf_mesh->getFaceCentroid(faces[i])[2] - subsurf_mesh->getFaceCentroid(faces[i + 1])[2];
      AMANZI_ASSERT(dz[0][cells[i]] > 0.);
    }
  }
  ATS::CLM::set_dz(dz);
  ATS::CLM::set_et_controls(1, 2, 0.1, 1.0, 0.1);
  ATS::CLM::setup_end();
  ATS::CLM::set_dz(dz);
}


void
SurfaceBalanceCLM::InitializePrimaryVariables_(const Tag& tag)
{
  // set as intialized the sources
  S_->GetW<CompositeVector>(surf_water_src_key_, tag, name_).PutScalar(0.);
  S_->GetRecordW(surf_water_src_key_, tag, name_).set_initialized();
  S_->GetW<CompositeVector>(ss_water_src_key_, tag, name_).PutScalar(0.);
  S_->GetRecordW(ss_water_src_key_, tag, name_).set_initialized();

  // set as intialized the diagnostics
  //  S_->GetRecordW(evap_flux_key_, tag, name_)->set_initialized();
  S_->GetRecordW(qE_lh_key_, tag, name_).set_initialized();
  S_->GetRecordW(qE_sh_key_, tag, name_).set_initialized();
  S_->GetRecordW(qE_lw_out_key_, tag, name_).set_initialized();
  S_->GetRecordW(qE_cond_key_, tag, name_).set_initialized();
  S_->GetRecordW(snow_swe_key_, tag, name_).set_initialized();
  S_->GetRecordW(can_wc_key_, tag, name_).set_initialized();
  S_->GetRecordW(surf_temp_key_, tag, name_).set_initialized();
  S_->GetRecordW(soil_temp_key_, tag, name_).set_initialized();
  S_->GetRecordW(can_temp_key_, tag, name_).set_initialized();
}


bool
SurfaceBalanceCLM::AdvanceStep(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  bool debug = false;
  Teuchos::RCP<VerboseObject> dcvo = Teuchos::null;
  int rank = mesh_->getComm()->MyPID();
  double dt = t_new - t_old;
  AMANZI_ASSERT(std::abs(dt - dt_) < 1.e-4);

  if (vo_->os_OK(Teuchos::VERB_LOW))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Advancing: t0 = " << S_->get_time(tag_current_)
               << " t1 = " << S_->get_time(tag_next_) << " h = " << dt << std::endl
               << "----------------------------------------------------------------" << std::endl;

  Tag tag = tag_current_;

  // Set the state
  S_->GetEvaluator(pres_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& pres =
    *S_->Get<CompositeVector>(pres_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(poro_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& poro =
    *S_->Get<CompositeVector>(poro_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(sl_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& sl =
    *S_->Get<CompositeVector>(sl_key_, tag).ViewComponent("cell", false);

  double patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  ATS::CLM::set_wc(poro, sl);
  ATS::CLM::set_tksat_from_porosity(poro);
  ATS::CLM::set_pressure(pres, patm);

  // set the forcing
  S_->GetEvaluator(met_sw_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& met_sw =
    *S_->Get<CompositeVector>(met_sw_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(met_lw_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& met_lw =
    *S_->Get<CompositeVector>(met_lw_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(met_air_temp_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& met_air_temp =
    *S_->Get<CompositeVector>(met_air_temp_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(met_vp_air_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& met_vp_air =
    *S_->Get<CompositeVector>(met_vp_air_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(met_wind_speed_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& met_wind_speed =
    *S_->Get<CompositeVector>(met_wind_speed_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(met_prain_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& met_prain =
    *S_->Get<CompositeVector>(met_prain_key_, tag).ViewComponent("cell", false);
  S_->GetEvaluator(met_psnow_key_, tag).Update(*S_, name_);
  const Epetra_MultiVector& met_psnow =
    *S_->Get<CompositeVector>(met_psnow_key_, tag).ViewComponent("cell", false);

  ATS::CLM::set_met_data(
    met_sw, met_lw, met_prain, met_psnow, met_air_temp, met_vp_air, met_wind_speed, patm);

  // set the start time, endtime
  ATS::CLM::advance_time(S_->get_cycle(tag), t_old, dt); // units in seconds

  // get diagnostics
  Epetra_MultiVector& qE_lh =
    *S_->GetW<CompositeVector>(qE_lh_key_, tag, name_).ViewComponent("cell", false);
  Epetra_MultiVector& qE_sh =
    *S_->GetW<CompositeVector>(qE_sh_key_, tag, name_).ViewComponent("cell", false);
  Epetra_MultiVector& qE_lw_out =
    *S_->GetW<CompositeVector>(qE_lw_out_key_, tag, name_).ViewComponent("cell", false);
  Epetra_MultiVector& qE_cond =
    *S_->GetW<CompositeVector>(qE_cond_key_, tag, name_).ViewComponent("cell", false);
  ATS::CLM::get_ground_energy_fluxes(qE_lh, qE_sh, qE_lw_out, qE_cond);

  Epetra_MultiVector& snow_depth =
    *S_->GetW<CompositeVector>(key_, tag, name_).ViewComponent("cell", false);
  Epetra_MultiVector& snow_swe =
    *S_->GetW<CompositeVector>(snow_swe_key_, tag, name_).ViewComponent("cell", false);
  Epetra_MultiVector& can_wc =
    *S_->GetW<CompositeVector>(can_wc_key_, tag, name_).ViewComponent("cell", false);
  Epetra_MultiVector& surf_temp =
    *S_->GetW<CompositeVector>(surf_temp_key_, tag, name_).ViewComponent("cell", false);
  Epetra_MultiVector& soil_temp =
    *S_->GetW<CompositeVector>(soil_temp_key_, tag, name_).ViewComponent("cell", false);
  Epetra_MultiVector& can_temp =
    *S_->GetW<CompositeVector>(can_temp_key_, tag, name_).ViewComponent("cell", false);
  ATS::CLM::get_diagnostics(snow_swe, snow_depth, can_wc, surf_temp, can_temp, soil_temp);
  changedEvaluatorPrimary(key_, tag, *S_);

  // get output
  Epetra_MultiVector& surf_water_src =
    *S_->GetW<CompositeVector>(surf_water_src_key_, tag, name_).ViewComponent("cell", false);
  Epetra_MultiVector& ss_water_src =
    *S_->GetW<CompositeVector>(ss_water_src_key_, tag, name_).ViewComponent("cell", false);
  ATS::CLM::get_total_mass_fluxes(surf_water_src, ss_water_src);
  changedEvaluatorPrimary(surf_water_src_key_, tag, *S_);
  changedEvaluatorPrimary(ss_water_src_key_, tag, *S_);

  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<std::string> vnames;
    std::vector<Teuchos::Ptr<const CompositeVector>> vecs;

    vnames.push_back("inc shortwave radiation [W/m^2]");
    vecs.push_back(S_->GetPtr<CompositeVector>(met_sw_key_, tag).ptr());
    vnames.push_back("inc longwave radiation [W/m^2]");
    vecs.push_back(S_->GetPtr<CompositeVector>(met_lw_key_, tag).ptr());
    vnames.push_back("inc latent heat [W/m^2]");
    vecs.push_back(S_->GetPtr<CompositeVector>(qE_lh_key_, tag).ptr());
    vnames.push_back("inc sensible heat [W/m^2]");
    vecs.push_back(S_->GetPtr<CompositeVector>(qE_sh_key_, tag).ptr());
    vnames.push_back("out longwave radiation [W/m^2]");
    vecs.push_back(S_->GetPtr<CompositeVector>(qE_lw_out_key_, tag).ptr());
    vnames.push_back("out conducted soil [W/m^2]");
    vecs.push_back(S_->GetPtr<CompositeVector>(qE_cond_key_, tag).ptr());

    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("surface water source [m/s]");
    vecs.push_back(S_->GetPtr<CompositeVector>(surf_water_src_key_, tag).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();

    vnames.clear();
    vecs.clear();

    vnames.push_back("snow depth [m]");
    vecs.push_back(S_->GetPtr<CompositeVector>(key_, tag).ptr());
    vnames.push_back("snow swe [m]");
    vecs.push_back(S_->GetPtr<CompositeVector>(snow_swe_key_, tag).ptr());
    vnames.push_back("canopy storage [m]");
    vecs.push_back(S_->GetPtr<CompositeVector>(can_wc_key_, tag).ptr());
    vnames.push_back("skin temperature [K]");
    vecs.push_back(S_->GetPtr<CompositeVector>(surf_temp_key_, tag).ptr());
    vnames.push_back("leaf temperature [K]");
    vecs.push_back(S_->GetPtr<CompositeVector>(can_temp_key_, tag).ptr());
    db_->WriteVectors(vnames, vecs, true);
    db_->WriteDivider();
  }
  return false;
}


} // namespace SurfaceBalance
} // namespace Amanzi
