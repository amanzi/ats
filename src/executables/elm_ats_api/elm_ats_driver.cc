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
createELM_ATSDriver(MPI_Fint *f_comm, const char *infile) {
  // -- create communicator & get process rank
  //auto comm = Amanzi::getDefaultComm();
  auto c_comm = MPI_Comm_f2c(*f_comm);
  auto comm = Amanzi::getComm(c_comm);
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
  return new ELM_ATSDriver(plist, comm);
}


ELM_ATSDriver::ELM_ATSDriver(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                             const Amanzi::Comm_ptr_type& comm)
  : Coordinator(plist, comm)
{
  // -- set default verbosity level to no output
  Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

  // domains names
  domain_subsurf_ = Amanzi::Keys::readDomain(*plist_, "domain");
  domain_surf_ = Amanzi::Keys::readDomainHint(*plist_, domain_subsurf_, "subsurface", "surface");

  // keys for fields exchanged with ELM
  infilt_key_ = Amanzi::Keys::readKey(*plist_, domain_subsurf_, "infiltration", "infiltration");
  pot_evap_key_ = Amanzi::Keys::readKey(*plist_, domain_subsurf_, "potential evaporation", "potential_evaporation");
  evap_key_ = Amanzi::Keys::readKey(*plist_, domain_subsurf_, "evaporation", "evaporation");
  trans_key_ = Amanzi::Keys::readKey(*plist_, domain_subsurf_, "transpiration", "transpiration");

  pres_key_ = Amanzi::Keys::readKey(*plist_, domain_subsurf_, "pressure", "pressure");
  pd_key_ = Amanzi::Keys::readKey(*plist_, domain_surf_, "ponded depth", "ponded_depth");
  satl_key_ = Amanzi::Keys::readKey(*plist_, domain_subsurf_, "saturation_liquid", "saturation_liquid");

  poro_key_ = Amanzi::Keys::readKey(*plist_, domain_subsurf_, "porosity", "base_porosity");
  elev_key_ = Amanzi::Keys::readKey(*plist_, domain_surf_, "elevation", "elevation");
  surf_cv_key_ = Amanzi::Keys::readKey(*plist_, domain_surf_, "surface cell volume", "cell_volume");

  // keys for fields used to convert ELM units to ATS units
  surf_mol_dens_key_ = Amanzi::Keys::readKey(*plist_, domain_surf_, "surface molar density", "molar_density_liquid");
  surf_mass_dens_key_ = Amanzi::Keys::readKey(*plist_, domain_surf_, "surface mass density", "mass_density_liquid");
  subsurf_mol_dens_key_ = Amanzi::Keys::readKey(*plist_, domain_subsurf_, "molar density", "molar_density_liquid");
  subsurf_mass_dens_key_ = Amanzi::Keys::readKey(*plist_, domain_subsurf_, "mass density", "mass_density_liquid");
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
  ncolumns_ = mesh_surf_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  AMANZI_ASSERT(ncolumns_ == mesh_subsurf_->num_columns(false));

  // -- get num cells per column - include consistency check later need to know
  //    if coupling zone is the entire subsurface mesh (as currently coded) or
  //    a portion of the total depth specified by # of cells into the
  //    subsurface
  auto& col_zero = mesh_subsurf_->cells_of_column(0);
  ncells_per_col_ = col_zero.size();

  // require primary variables (ELM --> ATS)
  requireAtNext(infilt_key_, Amanzi::Tags::NEXT, *S_, infilt_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  requireAtNext(pot_evap_key_, Amanzi::Tags::NEXT, *S_, pot_evap_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  requireAtNext(pot_trans_key_, Amanzi::Tags::NEXT, *S_, pot_trans_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);

  requireAtNext(poro_key_, Amanzi::Tags::NEXT, *S_, poro_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  requireAtNext(perm_key_, Amanzi::Tags::NEXT, *S_, perm_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);

  // require others (ATS --> ELM)
  requireAtNext(pd_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  requireAtNext(pres_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  requireAtNext(elev_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  requireAtNext(surf_cv_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);

  requireAtNext(surf_mol_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  requireAtNext(surf_mass_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_surf_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  requireAtNext(subsurf_mol_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  requireAtNext(subsurf_mass_dens_key_, Amanzi::Tags::NEXT, *S_)
    .SetMesh(mesh_subsurf_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  
  Coordinator::setup();
}


void ELM_ATSDriver::initialize(double t,
        double *p_atm,
        double *pressure)
{
  t0_ = t;

  // set initial conditions
  auto& pres = *S_->GetW<Amanzi::CompositeVector>(pres_key_, Amanzi::Tags::NEXT, pres_key_)
    .ViewComponent("cell", false);
  copySubsurf_(pres, pressure);
  // -- unclear whether to set initialized.  Need it to be false so that
  //    initializing faces is done via DeriveFaceValuesFromCellValues(), but it
  //    may never get set to true in that case?

  // set zero initial values for ELM exchange variables
  S_->GetW<Amanzi::CompositeVector>(infilt_key_, Amanzi::Tags::NEXT, infilt_key_)
    .PutScalar(0.);
  S_->GetRecordW(infilt_key_, Amanzi::Tags::NEXT, infilt_key_).set_initialized();

  S_->GetW<Amanzi::CompositeVector>(pot_evap_key_, Amanzi::Tags::NEXT, pot_evap_key_)
    .PutScalar(0.);
  S_->GetRecordW(pot_evap_key_, Amanzi::Tags::NEXT, pot_evap_key_).set_initialized();

  S_->GetW<Amanzi::CompositeVector>(pot_trans_key_, Amanzi::Tags::NEXT, pot_trans_key_)
    .PutScalar(0.);
  S_->GetRecordW(pot_trans_key_, Amanzi::Tags::NEXT, pot_trans_key_).set_initialized();

  // initialize fields, commit initial conditions
  Coordinator::initialize();
}


void ELM_ATSDriver::advance(double dt)
{
  double dt_subcycle = dt;
  double t_end = S_->get_time() + dt_subcycle;

  bool fail{false};
  while (S_->get_time() < t_end && dt_subcycle > 0.0) {
    S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt_subcycle);
    S_->advance_time(Amanzi::Tags::NEXT, dt_subcycle);

    // solve model for a single timestep
    fail = Coordinator::advance();

    if (fail) {
      // reset t_new
      S_->set_time(Amanzi::Tags::NEXT, S_->get_time(Amanzi::Tags::CURRENT));
    } else {
      S_->set_time(Amanzi::Tags::CURRENT, S_->get_time(Amanzi::Tags::NEXT));
      S_->advance_cycle();

      // make observations, vis, and checkpoints
      for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());
      visualize();
      checkpoint(); // checkpoint with the new dt
    }

    dt_subcycle = get_dt(fail);
  } // while

  if (fail) {
    // write one more vis for help debugging
    S_->advance_cycle(Amanzi::Tags::NEXT);
    visualize(true); // force vis

    // flush observations to make sure they are saved
    for (const auto& obs : observations_) obs->Flush();

    // dump a post_mortem checkpoint file for debugging
    checkpoint_->set_filebasename("post_mortem");
    checkpoint_->Write(*S_, Amanzi::Checkpoint::WriteType::POST_MORTEM);

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
    advance(dt);
  }
}

void ELM_ATSDriver::finalize()
{
  WriteStateStatistics(*S_, *vo_);
  report_memory();
  Teuchos::TimeMonitor::summarize(*vo_->os());
  Coordinator::finalize();
}

/*
assume that incoming data is in form
soil_infiltration(ncols)
soil_evaporation(ncols)
root_transpiration(ncells)
with ncells = ncells_per_column * ncols

assume surface can be indexed:
for (Amanzi::AmanziMesh::Entity_ID col=0; col!=ncolumns_; ++col)
{ ATS_surf_Epetra_MultiVector[0][col] = elm_surf_data[col];

  and assume subsurface can be indexed:
  auto& col_iter = mesh_sub_->cells_of_column(col);
  for (std::size_t i=0; i!=col_iter.size(); ++i)
  {  ATS_sub_Epetra_MultiVector[0][col_iter[i]] = elm_sub_data[col*ncells_per_col_+i]; }
}

assume evaporation and infiltration have same sign
*/
void
ELM_ATSDriver::set_potential_sources(double const *soil_infiltration,
        double const *soil_evaporation,
        double const *root_transpiration)
{
  // get densities to scale source fluxes
  S_->GetEvaluator(surf_mol_dens_key_, Amanzi::Tags::NEXT).Update(*S_, surf_mol_dens_key_);
  const Epetra_MultiVector& surf_mol_dens = *S_->Get<Amanzi::CompositeVector>(surf_mol_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  S_->GetEvaluator(surf_mass_dens_key_, Amanzi::Tags::NEXT).Update(*S_, surf_mass_dens_key_);
  const Epetra_MultiVector& surf_mass_dens = *S_->Get<Amanzi::CompositeVector>(surf_mass_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  S_->GetEvaluator(subsurf_mol_dens_key_, Amanzi::Tags::NEXT).Update(*S_, subsurf_mol_dens_key_);
  const Epetra_MultiVector& subsurf_mol_dens = *S_->Get<Amanzi::CompositeVector>(subsurf_mol_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  S_->GetEvaluator(subsurf_mass_dens_key_, Amanzi::Tags::NEXT).Update(*S_, subsurf_mass_dens_key_);
  const Epetra_MultiVector& subsurf_mass_dens = *S_->Get<Amanzi::CompositeVector>(subsurf_mass_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);

  // get sources
  Epetra_MultiVector& infilt = *S_->GetW<Amanzi::CompositeVector>(infilt_key_, Amanzi::Tags::NEXT, infilt_key_)
      .ViewComponent("cell", false);
  Epetra_MultiVector& pot_evap = *S_->GetW<Amanzi::CompositeVector>(pot_evap_key_, Amanzi::Tags::NEXT, pot_evap_key_)
      .ViewComponent("cell", false);
  Epetra_MultiVector& pot_trans = *S_->GetW<Amanzi::CompositeVector>(pot_trans_key_, Amanzi::Tags::NEXT, pot_trans_key_)
      .ViewComponent("cell", false);

  // scale evaporation and infiltration and add to surface source
  // assume evaporation and infiltration have same sign?
  // negative out of subsurface and positive into subsurface?
  // or positive evap is out and positive infil is in?
  // for now, assume same sign
  for (Amanzi::AmanziMesh::Entity_ID col=0; col!=ncolumns_; ++col) {
    double surf_mol_h20_kg = surf_mol_dens[0][col] / surf_mass_dens[0][col];
    infilt[0][col] = soil_infiltration[col] * surf_mol_h20_kg;
    pot_evap[0][col] = soil_evaporation[col] * surf_mol_h20_kg;

    // scale root_transpiration and add to subsurface source
    auto& col_iter = mesh_subsurf_->cells_of_column(col);
    for (std::size_t i=0; i!=col_iter.size(); ++i) {
      double subsurf_mol_h20_kg = subsurf_mol_dens[0][col_iter[i]] / subsurf_mass_dens[0][col_iter[i]];
      pot_trans[0][col_iter[i]] = root_transpiration[col*ncells_per_col_+i] * subsurf_mol_h20_kg;
    }
  }

  // mark sources as changed
  changedEvaluatorPrimary(infilt_key_, Amanzi::Tags::NEXT, *S_);
  changedEvaluatorPrimary(pot_evap_key_, Amanzi::Tags::NEXT, *S_);
  changedEvaluatorPrimary(pot_trans_key_, Amanzi::Tags::NEXT, *S_);
}


void
ELM_ATSDriver::get_actual_sources(double *soil_infiltration,
        double *soil_evaporation,
        double *root_transpiration)
{
  // get densities to scale source fluxes
  S_->GetEvaluator(surf_mol_dens_key_, Amanzi::Tags::NEXT).Update(*S_, surf_mol_dens_key_);
  const Epetra_MultiVector& surf_mol_dens = *S_->Get<Amanzi::CompositeVector>(surf_mol_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  S_->GetEvaluator(surf_mass_dens_key_, Amanzi::Tags::NEXT).Update(*S_, surf_mass_dens_key_);
  const Epetra_MultiVector& surf_mass_dens = *S_->Get<Amanzi::CompositeVector>(surf_mass_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  S_->GetEvaluator(subsurf_mol_dens_key_, Amanzi::Tags::NEXT).Update(*S_, subsurf_mol_dens_key_);
  const Epetra_MultiVector& subsurf_mol_dens = *S_->Get<Amanzi::CompositeVector>(subsurf_mol_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  S_->GetEvaluator(subsurf_mass_dens_key_, Amanzi::Tags::NEXT).Update(*S_, subsurf_mass_dens_key_);
  const Epetra_MultiVector& subsurf_mass_dens = *S_->Get<Amanzi::CompositeVector>(subsurf_mass_dens_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);

  // get sources
  const Epetra_MultiVector& infilt = *S_->GetW<Amanzi::CompositeVector>(infilt_key_, Amanzi::Tags::NEXT, infilt_key_)
      .ViewComponent("cell", false);
  const Epetra_MultiVector& pot_evap = *S_->GetW<Amanzi::CompositeVector>(pot_evap_key_, Amanzi::Tags::NEXT, pot_evap_key_)
      .ViewComponent("cell", false);
  const Epetra_MultiVector& pot_trans = *S_->GetW<Amanzi::CompositeVector>(pot_trans_key_, Amanzi::Tags::NEXT, pot_trans_key_)
      .ViewComponent("cell", false);

  // scale evaporation and infiltration and add to surface source
  // assume evaporation and infiltration have same sign?
  // negative out of subsurface and positive into subsurface?
  // or positive evap is out and positive infil is in?
  // for now, assume same sign
  for (Amanzi::AmanziMesh::Entity_ID col=0; col!=ncolumns_; ++col) {
    double surf_mol_h20_kg = surf_mol_dens[0][col] / surf_mass_dens[0][col];
    soil_infiltration[col] = infilt[0][col] / surf_mol_h20_kg;
    soil_evaporation[col] = pot_evap[0][col] / surf_mol_h20_kg;

    // scale root_transpiration and add to subsurface source
    auto& col_iter = mesh_subsurf_->cells_of_column(col);
    for (std::size_t i=0; i!=col_iter.size(); ++i) {
      double subsurf_mol_h20_kg = subsurf_mol_dens[0][col_iter[i]] / subsurf_mass_dens[0][col_iter[i]];
      root_transpiration[col*ncells_per_col_+i] = pot_trans[0][col_iter[i]] / subsurf_mol_h20_kg;
    }
  }
}


//
// soil_potential and sat_ice are currently ignored/set to zero
//
void
ELM_ATSDriver::get_waterstate(double *ponded_depth,
        double *soil_pressure,
        double *soil_potential,
        double *sat_liq,
        double *sat_ice)
{
  S_->GetEvaluator(pd_key_, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  const auto& pd = *S_->Get<Amanzi::CompositeVector>(pd_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  S_->GetEvaluator(pres_key_, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  const auto& pres = *S_->Get<Amanzi::CompositeVector>(pres_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);
  S_->GetEvaluator(satl_key_, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  const auto& sat = *S_->Get<Amanzi::CompositeVector>(satl_key_, Amanzi::Tags::NEXT)
      .ViewComponent("cell", false);

  for (Amanzi::AmanziMesh::Entity_ID col=0; col!=ncolumns_; ++col) {
    ponded_depth[col] = pd[0][col];

    auto& col_iter = mesh_subsurf_->cells_of_column(col);
    for (std::size_t i=0; i!=col_iter.size(); ++i) {
      soil_pressure[col*ncells_per_col_+i] = pres[0][col_iter[i]];
      soil_potential[col*ncells_per_col_+i] = std::max(0.0, 101325. - pres[0][col_iter[i]]);
      sat_liq[col*ncells_per_col_+i] = sat[0][col_iter[i]];
      sat_ice[col*ncells_per_col_+i] = 0.;
    }
  }
}


//
// dz & depth are currently ignored -- presumed identical between ATS & ELM
//
void ELM_ATSDriver::get_mesh_info(int& ncols_local,
        int& ncols_global,
        int& ncells_per_col,
        double* lat,
        double* lon,
        double* elev,
        double* surface_area,
        double* dz,
        double* depth)
{
  ncols_local = ncolumns_;
  ncols_global = mesh_surf_->cell_map(Amanzi::AmanziMesh::Entity_kind::CELL).NumGlobalElements();
  ncells_per_col = ncells_per_col_;

  S_->GetEvaluator(elev_key_, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  copySurf_(elev, *S_->Get<Amanzi::CompositeVector>(elev_key_, Amanzi::Tags::NEXT)
            .ViewComponent("cell", false));

  S_->GetEvaluator(surf_cv_key_, Amanzi::Tags::NEXT).Update(*S_, "ELM");
  copySurf_(surface_area, *S_->Get<Amanzi::CompositeVector>(surf_cv_key_, Amanzi::Tags::NEXT)
            .ViewComponent("cell", false));

  // hard-coded Toledo OH for now...
  copySurf_(lat, 41.65);
  copySurf_(lon, -83.54);
}


} // namespace ATS
