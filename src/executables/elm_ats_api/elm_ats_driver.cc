#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include "errors.hh"
#include "dbc.hh"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "VerboseObject.hh"

// registration files
#include "state_evaluators_registration.hh"
#include "ats_relations_registration.hh"
#include "ats_transport_registration.hh"
#include "ats_energy_pks_registration.hh"
#include "ats_energy_relations_registration.hh"
#include "ats_flow_pks_registration.hh"
#include "ats_flow_relations_registration.hh"
#include "ats_deformation_registration.hh"
#include "ats_bgc_registration.hh"
#include "ats_surface_balance_registration.hh"
#include "ats_mpc_registration.hh"
//#include "ats_sediment_transport_registration.hh"
#include "mdm_transport_registration.hh"
#include "multiscale_transport_registration.hh"
#ifdef ALQUIMIA_ENABLED
#include "pks_chemistry_registration.hh"
#endif

// include fenv if it exists
#include "boost/version.hpp"
#if (BOOST_VERSION / 100 % 1000 >= 46)
#include "boost/config.hpp"
#ifndef BOOST_NO_FENV_H
#ifdef _GNU_SOURCE
#define AMANZI_USE_FENV
#include "boost/detail/fenv.hpp"
#endif
#endif
#endif
#include "boost/filesystem.hpp"

#include "AmanziComm.hh"
#include "CompositeVector.hh"
#include "IO.hh"
#include "UnstructuredObservations.hh"
#include "pk_helpers.hh"
#include "elm_ats_driver.hh"

namespace ATS {

ELM_ATSDriver*
createELM_ATSDriver(MPI_Fint *f_comm, const char *infile, int npfts) {
  // -- create communicator & get process rank
  //auto comm = getDefaultComm();
  auto c_comm = MPI_Comm_f2c(*f_comm);
  auto comm = getComm(c_comm);
  auto rank = comm->MyPID();

  // convert input file to std::string for easier handling
  // infile must be null-terminated
  std::string input_filename(infile);

  // check validity of input file name
  if (input_filename.empty()) {
    if (rank == 0)
      std::cerr << "ERROR: no input file provided" << std::endl;
  } else if (!boost::filesystem::exists(input_filename)) {
    if (rank == 0)
      std::cerr << "ERROR: input file \"" << input_filename << "\" does not exist." << std::endl;
  }

  // -- parse input file
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(input_filename);
  return new ELM_ATSDriver(plist, comm, npfts);
}


ELM_ATSDriver::ELM_ATSDriver(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                             const Comm_ptr_type& comm,
                             int npfts)
  : Coordinator(plist, comm),
    npfts_(npfts),
    ncolumns_(-1),
    ncells_per_col_(-1)
{
  // -- set default verbosity level to no output
  VerboseObject::global_default_level = Teuchos::VERB_NONE;

  // domains names
  domain_subsurf_ = Keys::readDomain(*plist_, "domain");
  domain_surf_ = Keys::readDomainHint(*plist_, domain_subsurf_, "subsurface", "surface");

  // keys for fields in mesh info
  // lat_key_ = Keys::readKey(*plist_, domain_surf_, "latitude", "latitude");
  // lon_key_ = Keys::readKey(*plist_, domain_surf_, "longitude", "longitude");
  elev_key_ = Keys::readKey(*plist_, domain_surf_, "elevation", "elevation");

  // keys for initialize
  // pres_atm_key_ = Keys::readKey(*plist_, domain_surf_, "atmospheric_pressure", "atmospheric_pressure");
  pres_key_ = Keys::readKey(*plist_, domain_subsurf_, "pressure", "pressure");

  // soil parameters
  base_poro_key_ = Keys::readKey(*plist_, domain_subsurf_, "base porosity", "base_porosity");
  perm_key_ = Keys::readKey(*plist_, domain_subsurf_, "permeability", "permeability");
  // ch_b_key_ = Keys::readKey(*plist_, domain_subsurf_, "Clapp and Hornberger b", "clapp_horn_b");
  // ch_smpsat_key_ = Keys::readKey(*plist_, domain_subsurf_, "Clapp and Hornberger soil mafic potential at saturation", "clapp_horn_smpsat");
  // ch_sr_key_ = Keys::readKey(*plist_, domain_subsurf_, "Clapp and Hornberger residual saturation", "clapp_horn_sr");

  // soil properties
  // poro_key_ = Keys::readKey(*plist_, domain_subsurf_, "porosity", "porosity");
  root_frac_key_ = Keys::readKey(*plist_, domain_subsurf_, "rooting fraction", "rooting_fraction");

  // potential sources
  pot_infilt_key_ = Keys::readKey(*plist_, domain_surf_, "potential infiltration", "potential_infiltration");
  pot_evap_key_ = Keys::readKey(*plist_, domain_surf_, "potential evaporation", "potential_evaporation");
  pot_trans_key_ = Keys::readKey(*plist_, domain_surf_, "potential transpiration", "potential_transpiration");

  // water state
  pd_key_ = Keys::readKey(*plist_, domain_surf_, "ponded depth", "ponded_depth");
  pc_key_ = Keys::readKey(*plist_, domain_surf_, "capillary pressure, air over liquid", "capillary_pressure_air_liq");
  sat_liq_key_ = Keys::readKey(*plist_, domain_subsurf_, "saturation liquid", "saturation_liquid");
  sat_ice_key_ = Keys::readKey(*plist_, domain_subsurf_, "saturation ice", "saturation_ice");

  // water fluxes
  infilt_key_ = Keys::readKey(*plist_, domain_surf_, "infiltration", "surface_subsurface_flux");
  evap_key_ = Keys::readKey(*plist_, domain_surf_, "evaporation", "evaporation");
  trans_key_ = Keys::readKey(*plist_, domain_subsurf_, "transpiration", "transpiration");

  // keys for fields used to convert ELM units to ATS units
  surf_mol_dens_key_ = Keys::readKey(*plist_, domain_surf_, "surface molar density", "molar_density_liquid");
  surf_mass_dens_key_ = Keys::readKey(*plist_, domain_surf_, "surface mass density", "mass_density_liquid");
  subsurf_mol_dens_key_ = Keys::readKey(*plist_, domain_subsurf_, "molar density", "molar_density_liquid");
  subsurf_mass_dens_key_ = Keys::readKey(*plist_, domain_subsurf_, "mass density", "mass_density_liquid");

  // cell vol keys
  surf_cv_key_ = Keys::getKey(domain_surf_, "cell_volume");
  cv_key_ = Keys::getKey(domain_subsurf_, "cell_volume");
}


void
ELM_ATSDriver::setup()
{
  // meshing
  mesh_subsurf_ = S_->GetMesh(domain_subsurf_);
  mesh_surf_ = S_->GetMesh(domain_surf_);

  // -- build columns to allow indexing by column
  mesh_subsurf_->build_columns();

  // -- check that number of surface cells = number of columns
  ncolumns_ = mesh_surf_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  AMANZI_ASSERT(ncolumns_ == mesh_subsurf_->num_columns(false));

  // -- get num cells per column - include consistency check later need to know
  //    if coupling zone is the entire subsurface mesh (as currently coded) or
  //    a portion of the total depth specified by # of cells into the
  //    subsurface
  auto& col_zero = mesh_subsurf_->cells_of_column(0);
  ncells_per_col_ = col_zero.size();
  for (int col=0; col!=ncolumns_; ++col)
    AMANZI_ASSERT(mesh_subsurf_->cells_of_column(col).size() == ncells_per_col_);

  // require primary variables (ELM --> ATS)
  requireAtNext(base_poro_key_, Tags::NEXT, *S_, base_poro_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(perm_key_, Tags::NEXT, *S_, perm_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(root_frac_key_, Tags::NEXT, *S_, perm_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(pot_infilt_key_, Tags::NEXT, *S_, pot_infilt_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(pot_evap_key_, Tags::NEXT, *S_, pot_evap_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(pot_trans_key_, Tags::NEXT, *S_, pot_trans_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", AmanziMesh::CELL, 1);

  // require others (ATS --> ELM)
  // requireAtNext(lat_key_, Tags::NEXT, *S_)
  //   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  // requireAtNext(lon_key_, Tags::NEXT, *S_)
  //   .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(elev_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(surf_cv_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(pd_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(pres_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(pc_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(sat_liq_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(sat_ice_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(infilt_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(evap_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(trans_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(surf_mol_dens_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(surf_mass_dens_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(subsurf_mol_dens_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  requireAtNext(subsurf_mass_dens_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);

  requireAtNext(cv_key_, Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", AmanziMesh::CELL, 1);
  Coordinator::setup();
}



//
// dz & depth are currently ignored -- presumed identical between ATS & ELM
//
void ELM_ATSDriver::get_mesh_info(int& ncols_local,
        int& ncols_global,
        double * const lat,
        double * const lon,
        double * const elev,
        double * const surface_area,
        int * const pft,
        int& nlevgrnd,
        double * const depth)
{
  ncols_local = ncolumns_;
  ncols_global = mesh_surf_->cell_map(AmanziMesh::Entity_kind::CELL).NumGlobalElements();

  copyFromSurf_(elev, elev_key_);
  copyFromSurf_(surface_area, surf_cv_key_);
  // copyFromSurf_(lat, lat_key_);
  // copyFromSurf_(lon, lon_key_);

  // NOTE: figure out how to map from ATS LC types to ELM PFT... --etc
  // hard-coded veg type for now...
  for (int i=0; i!=ncolumns_; ++i) pft[i] = 1;

  nlevgrnd = ncells_per_col_;
  const auto& cells_in_col = mesh_subsurf_->cells_of_column(0);
  const auto& fc = mesh_subsurf_->face_centroid(mesh_subsurf_->faces_of_column(0)[0]);
  for (int i=0; i!=ncells_per_col_; ++i) {
    depth[i] = fc[2] - mesh_subsurf_->cell_centroid(cells_in_col[i])[2];
  }

  // hard-coded Toledo OH for now...
  for (int i=0; i!=ncolumns_; ++i) {
    lat[i] = 41.65;
    lon[i] = -83.54;
  }
}


void ELM_ATSDriver::initialize(double t,
        double const * const p_atm,
        double const * const pressure)
{
  t0_ = t;

  // set initial conditions on pressure:
  // -- set subsurface cells
  copyToSub_(pressure, pres_key_);
  // -- set subsurface faces
  DeriveFaceValuesFromCellValues(S_->GetW<CompositeVector>(pres_key_, Tags::NEXT, pres_key_));
  // -- tag as initialized
  S_->GetRecordW(pres_key_, Tags::NEXT, pres_key_).set_initialized();

  // set as zero and tag as initialized
  initZero_(pot_infilt_key_);
  initZero_(pot_trans_key_);
  initZero_(pot_evap_key_);

  // initialize ATS data, commit initial conditions
  Coordinator::initialize();
}


void ELM_ATSDriver::advance(double dt, bool do_chkp, bool do_vis)
{
  double dt_subcycle = dt;
  double t_end = S_->get_time() + dt_subcycle;

  bool fail{false};
  while (S_->get_time() < t_end && dt_subcycle > 0.0) {
    S_->Assign<double>("dt", Tags::DEFAULT, "dt", dt_subcycle);
    S_->advance_time(Tags::NEXT, dt_subcycle);

    // solve model for a single timestep
    fail = Coordinator::advance();

    if (fail) {
      // reset t_new
      S_->set_time(Tags::NEXT, S_->get_time(Tags::CURRENT));
    } else {
      S_->set_time(Tags::CURRENT, S_->get_time(Tags::NEXT));
      S_->advance_cycle();

      // make observations, vis, and checkpoints
      for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());

      // vis/checkpoint if EITHER ATS or ELM request it
      if (do_vis && !visualize()) visualize(true);
      if (do_chkp && !checkpoint()) checkpoint(true);
    }

    dt_subcycle = get_dt(fail);
  } // while

  if (fail) {
    // write one more vis for help debugging
    S_->advance_cycle(Tags::NEXT);
    visualize(true); // force vis

    // flush observations to make sure they are saved
    for (const auto& obs : observations_) obs->Flush();

    // dump a post_mortem checkpoint file for debugging
    checkpoint_->set_filebasename("post_mortem");
    checkpoint_->Write(*S_, Checkpoint::WriteType::POST_MORTEM);

    Errors::Message msg("ELM_ATSDriver: advance(dt) failed.");
    Exceptions::amanzi_throw(msg);
  }
} // advance()

// simulates external timeloop with dt coming from calling model
void ELM_ATSDriver::advance_test()
{
  while (S_->get_time() < t1_) {
    // use dt from ATS for testing
    double dt = Coordinator::get_dt(false);
    // call main method
    advance(dt, false, false);
  }
}

void ELM_ATSDriver::finalize()
{
  WriteStateStatistics(*S_, *vo_);
  report_memory();
  Teuchos::TimeMonitor::summarize(*vo_->os());
  Coordinator::finalize();
}


void ELM_ATSDriver::set_soil_hydrologic_parameters(double const * const base_porosity,
        double const * const hydraulic_conductivity,
        double const * const clapp_horn_b,
        double const * const clapp_horn_smpsat,
        double const * const clapp_horn_sr)
{
  copyToSub_(base_porosity, base_poro_key_);

  // convert Ksat to perm via rho * g / visc using default rho and visc values.
  copyToSub_(hydraulic_conductivity, perm_key_);
  auto& perm = *S_->GetW<CompositeVector>(perm_key_, Tags::NEXT, perm_key_)
    .ViewComponent("cell", false);
  double factor = 8.9e-4 / (1000 * -S_->Get<AmanziGeometry::Point>("gravity", Tags::NEXT)[2]);
  perm.Scale(factor);

  // copyToSub_(clapp_horn_b, ch_b_key_);
  // copyToSub_(clapp_horn_smpsat, ch_smpsat_key_);
  // copyToSub_(clapp_horn_sr, ch_sr_key_);
}

void ELM_ATSDriver::set_veg_parameters(double const * const mafic_potential_full_turgor,
        double const * const mafic_potential_wilt_point)
{
  // pass for now! FIXME --etc
}


void ELM_ATSDriver::set_soil_hydrologic_properties(double const * const effective_porosity)
{
  // this isn't well defined -- pass for now --etc
}


void ELM_ATSDriver::set_veg_properties(double const * const rooting_fraction)
{
  copyToSub_(rooting_fraction, root_frac_key_);
}


void
ELM_ATSDriver::set_potential_sources(double const * const infiltration,
        double const * const evaporation,
        double const * const transpiration)
{
  // note, this assumes that all are in units of m/s
  copyToSurf_(infiltration, pot_infilt_key_);
  copyToSurf_(evaporation, pot_evap_key_);
  copyToSurf_(transpiration, pot_trans_key_);
}


//
// soil_potential and sat_ice are currently ignored/set to zero
//
void
ELM_ATSDriver::get_waterstate(double * const ponded_depth,
        double * const soil_pressure,
        double * const soil_potential,
        double * const sat_liq,
        double * const sat_ice)
{
  copyFromSurf_(ponded_depth, pd_key_);
  copyFromSub_(soil_pressure, pres_key_);
  copyFromSub_(soil_potential, pc_key_);
  copyFromSub_(sat_liq, sat_liq_key_);
  // copyFromSub_(sat_ice, sat_ice_key_);
}


void
ELM_ATSDriver::get_water_fluxes(double * const infiltration,
        double * const evaporation,
        double * const transpiration,
        double * const net_subsurface_fluxes,
        double * const net_runon)
{
  copyFromSurf_(infiltration, infilt_key_); // confirm what this is used for --
                                            // pot_infilt is assumed to be net
                                            // surface fluxes not include evap,
                                            // while this is actual
                                            // infiltration --ETC
  copyFromSurf_(evaporation, evap_key_);

  // convert trans from mol/m^3/s to... ? --ETC
  copyFromSurf_(transpiration, trans_key_);

  // unclear how to implement net_subsurface_fluxes and net_runon... --ETC
}


void ELM_ATSDriver::initZero_(const Key& key)
{
  auto& vec = S_->GetW<CompositeVector>(key, Tags::NEXT, key);
  vec.PutScalar(0.);
  S_->GetRecordW(key, Tags::NEXT, key).set_initialized();
}


void ELM_ATSDriver::copyToSurf_(double const * const in, const Key& key)
{
  // surf maps directly into columns
  auto& vec = *S_->GetW<CompositeVector>(key, Tags::NEXT, key)
    .ViewComponent("cell", false);
  AMANZI_ASSERT(vec.MyLength() == ncolumns_);

  for (int i=0; i!=ncolumns_; ++i) vec[0][i] = in[i];

  changedEvaluatorPrimary(key, Tags::NEXT, *S_);
}


//
// ELM data is defined as var(col,lev), meaning that, due to Fortran ordering,
// the column is fasted varying, not the grid cell.
//
void ELM_ATSDriver::copyToSub_(double const * const in, const Key& key)
{
  auto& vec = *S_->GetW<CompositeVector>(key, Tags::NEXT, key)
    .ViewComponent("cell", false);

  for (int i=0; i!=ncolumns_; ++i) {
    const auto& cells_of_col = mesh_subsurf_->cells_of_column(i);
    for (int j=0; j!=ncells_per_col_; ++j) {
      vec[0][cells_of_col[j]] = in[j * ncolumns_ + i];
    }
  }

  changedEvaluatorPrimary(key, Tags::NEXT, *S_);
}


void ELM_ATSDriver::copyFromSurf_(double * const out, const Key& key) const
{
  S_->GetEvaluator(key, Tags::NEXT).Update(*S_, "ELM");
  const auto& vec = *S_->Get<CompositeVector>(key, Tags::NEXT)
    .ViewComponent("cell", false);

  for (int i=0; i!=ncolumns_; ++i) out[i] = vec[0][i];
}


//
// ELM data is defined as var(col,lev), meaning that, due to Fortran ordering,
// the column is fasted varying, not the grid cell.
//
void ELM_ATSDriver::copyFromSub_(double * const out, const Key& key) const
{
  S_->GetEvaluator(key, Tags::NEXT).Update(*S_, "ELM");
  const auto& vec = *S_->Get<CompositeVector>(key, Tags::NEXT)
    .ViewComponent("cell", false);

  for (int i=0; i!=ncolumns_; ++i) {
    const auto& cells_of_col = mesh_subsurf_->cells_of_column(i);
    for (int j=0; j!=ncells_per_col_; ++j) {
      out[j * ncolumns_ + i] = vec[0][cells_of_col[j]];
    }
  }
}


} // namespace ATS
