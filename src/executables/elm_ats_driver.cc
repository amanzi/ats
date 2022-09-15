
#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "VerboseObject_objs.hh"
#include "VerboseObject.hh"

#include "ats_version.hh"
#include "amanzi_version.hh"
#include "tpl_versions.h"

#include "dbc.hh"
#include "errors.hh"
#include "ats_driver.hh"

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
#include "AmanziTypes.hh"
#include "GeometricModel.hh"
#include "State.hh"

#include "pk_helpers.hh"
#include "ats_mesh_factory.hh"
#include "elm_ats_coordinator.hh"
#include "elm_ats_driver.hh"

namespace ATS {


  ELM_ATSDriver::ELM_ATSDriver()
  : Coordinator() {}

void
ELM_ATSDriver::setup(MPI_Fint *f_comm, const char *infile)
{
  // -- create communicator & get process rank
  //auto comm = Amanzi::getDefaultComm();
  auto c_comm = MPI_Comm_f2c(*f_comm);
  auto comm = setComm(c_comm);
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

  // instantiates much of Coordinator object
  // it is very import that this is called
  Coordinator:coordinator_init(*plist, comm)

  // -- set default verbosity level to no output
  Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

  // domains
  domain_sub_ = plist->get<std::string>("domain name", "domain");
  domain_srf_ = Amanzi::Keys::readDomainHint(*plist, domain_sub_, "subsurface", "surface");

  // keys for fields exchanged with ELM
  sub_src_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "subsurface source", "source_sink");
  srf_src_key_ = Amanzi::Keys::readKey(*plist, domain_srf_, "surface source", "source_sink");
  pres_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "pressure", "pressure");
  pd_key_ = Amanzi::Keys::readKey(*plist, domain_srf_, "ponded depth", "ponded_depth");
  satl_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "saturation_liquid", "saturation_liquid");
  por_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "porosity", "porosity");
  elev_key_ = Amanzi::Keys::readKey(*plist, domain_srf_, "elevation", "elevation");

  // keys for fields used to convert ELM units to ATS units
  srf_mol_dens_key_ = Amanzi::Keys::readKey(*plist, domain_srf_, "surface molar density", "molar_density_liquid");
  srf_mass_dens_key_ = Amanzi::Keys::readKey(*plist, domain_srf_, "surface mass density", "mass_density_liquid");
  sub_mol_dens_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "molar density", "molar_density_liquid");
  sub_mass_dens_key_ = Amanzi::Keys::readKey(*plist, domain_sub_, "mass density", "mass_density_liquid");

  // assume for now that mesh info has been communicated
  mesh_subsurf_ = S_->GetMesh(domain_sub_);
  mesh_surf_ = S_->GetMesh(domain_srf_);

  // build columns to allow indexing by column
  mesh_subsurf_->build_columns();

  // check that number of surface cells = number of columns
  ncolumns_ = mesh_surf_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  AMANZI_ASSERT(ncolumns_ == mesh_subsurf_->num_columns(false));

  // get num cells per column - include consistency check later
  // need to know if coupling zone is the entire subsurface mesh (as currently coded)
  // or a portion of the total depth specified by # of cells into the subsurface
  auto& col_zero = mesh_subsurf_->cells_of_column(0);
  ncol_cells_ = col_zero.size();

  // require primary variables
  // -- subsurface water source
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(sub_src_key_, Amanzi::Tags::NEXT,  sub_src_key_)
    .SetMesh(mesh_subsurf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  RequireEvaluatorPrimary(sub_src_key_, Amanzi::Tags::NEXT, *S_);
  // -- surface water source-sink
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(srf_src_key_, Amanzi::Tags::NEXT,  srf_src_key_)
    .SetMesh(mesh_surf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  RequireEvaluatorPrimary(srf_src_key_, Amanzi::Tags::NEXT, *S_);
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>("dz", Amanzi::Tags::NEXT,  "dz")
    .SetMesh(mesh_subsurf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  RequireEvaluatorPrimary("dz", Amanzi::Tags::NEXT, *S_);
  S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>("depth", Amanzi::Tags::NEXT,  "depth")
    .SetMesh(mesh_subsurf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  RequireEvaluatorPrimary("depth", Amanzi::Tags::NEXT, *S_);

  // commented out for now while building and testing- setting as Primary (and initializing as 0.0)
  // causes problems when running existing unmodified ATS xml files
  // -- porosity
  //S_->Require<Amanzi::CompositeVector,Amanzi::CompositeVectorSpace>(por_key_, Amanzi::Tags::NEXT,  por_key_)
  //  .SetMesh(mesh_subsurf_)->SetComponent("cell", Amanzi::AmanziMesh::CELL, 1);
  //RequireEvaluatorPrimary(por_key_, Amanzi::Tags::NEXT, *S_);

  // setup time tags and call setup() for PKs and state
  //Teuchos::TimeMonitor monitor(*setup_timer_);
  Coordinator::setup();
}


void ELM_ATSDriver::initialize()
{
  // initialize fields, commit initial conditions
  Coordinator::initialize();

  // set initial values for ELM exchange variables
  S_->GetW<Amanzi::CompositeVector>(sub_src_key_, Amanzi::Tags::NEXT, sub_src_key_).PutScalar(0.);
  S_->GetRecordW(sub_src_key_, Amanzi::Tags::NEXT, sub_src_key_).set_initialized();
  S_->GetW<Amanzi::CompositeVector>(srf_src_key_, Amanzi::Tags::NEXT, srf_src_key_).PutScalar(0.);
  S_->GetRecordW(srf_src_key_, Amanzi::Tags::NEXT, srf_src_key_).set_initialized();
  S_->GetW<Amanzi::CompositeVector>("dz", Amanzi::Tags::NEXT, "dz").PutScalar(0.);
  S_->GetRecordW("dz", Amanzi::Tags::NEXT, "dz").set_initialized();
  S_->GetW<Amanzi::CompositeVector>("depth", Amanzi::Tags::NEXT, "depth").PutScalar(0.);
  S_->GetRecordW("depth", Amanzi::Tags::NEXT, "depth").set_initialized();
  //S_->GetW<Amanzi::CompositeVector>(por_key_, Amanzi::Tags::NEXT, por_key_).PutScalar(0.);
  //S_->GetRecordW(por_key_, Amanzi::Tags::NEXT, por_key_).set_initialized();

  // visualization at IC
  // for testing
  visualize();
  checkpoint();

  // Make sure times are set up correctly
  AMANZI_ASSERT(std::abs(S_->get_time(Amanzi::Tags::NEXT)
                         - S_->get_time(Amanzi::Tags::CURRENT)) < 1.e-4);
}


void ELM_ATSDriver::advance(double *dt)
{
  double dt_subcycle = *dt;

  // should adapt this to use Coordinator::advance() to solve single steps
  //auto fail = Coordinator::advance(*dt);
  // start and end times for timestep
  double t_end = S_->get_time() + dt_subcycle;

  bool fail = false;
  while (S_->get_time() < t_end && dt_subcycle > 0.0) {

    // run model for a duration of dt
    // order is important
    // advance NEXT time tag
    S_->advance_time(Amanzi::Tags::NEXT, dt_subcycle);
    // get time from tags
    double t_old = S_->get_time(Amanzi::Tags::CURRENT);
    double t_new = S_->get_time(Amanzi::Tags::NEXT);
    // check that dt and time tags align
    AMANZI_ASSERT(std::abs((t_new - t_old) - dt_subcycle) < 1.e-4);

    // advance pks
    fail = pk_->AdvanceStep(t_old, t_new, false);
    if (!fail) fail |= !pk_->ValidStep();

    if (fail) {

      // Failed the timestep.
      // Potentially write out failed timestep for debugging
      for (const auto& vis : failed_visualization_) WriteVis(*vis, *S_);

      // set as failed and revert timestamps
      pk_->FailStep(t_old, t_new, Amanzi::Tags::NEXT);
      // set NEXT = CURRENT to reset t_new
      S_->set_time(Amanzi::Tags::NEXT, S_->get_time(Amanzi::Tags::CURRENT));

    } else {

      // commit the state, copying NEXT --> CURRENT
      pk_->CommitStep(t_old, t_new, Amanzi::Tags::NEXT);

      // set CURRENT = NEXT
      S_->set_time(Amanzi::Tags::CURRENT, S_->get_time(Amanzi::Tags::NEXT));
      // state advance_cycle - necessary??
      S_->advance_cycle();

      // make observations, vis, and checkpoints
      for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());
      visualize();
      checkpoint(); // checkpoint with the new dt
    }

    // get new dt and assign to State
    dt_subcycle = get_dt(fail);
    S_->Assign<double>("dt", Amanzi::Tags::DEFAULT, "dt", dt_subcycle);
  } // end while

  if (fail) {
    Errors::Message msg("ELM_ATSDriver: advance(dt) failed.");
    Exceptions::amanzi_throw(msg);
  }

  // update ATS->ELM data if necessary
  S_->GetEvaluator(pres_key_, Amanzi::Tags::NEXT).Update(*S_, pres_key_);
  const Epetra_MultiVector& pres = *S_->Get<Amanzi::
    CompositeVector>(pres_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);
  S_->GetEvaluator(satl_key_, Amanzi::Tags::NEXT).Update(*S_, satl_key_);
  const Epetra_MultiVector& satl = *S_->Get<Amanzi::
    CompositeVector>(satl_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);
  //S_->GetEvaluator(por_key_, Amanzi::Tags::NEXT).Update(*S_, por_key_);
  //const Epetra_MultiVector& poro = *S_->Get<Amanzi::CompositeVector>(por_key_, Amanzi::Tags::NEXT).ViewComponent("cell", false);
}

// simulates external timeloop with dt coming from calling model
void ELM_ATSDriver::advance_test()
{
  while (S_->get_time() < elm_coordinator_->get_end_time()) {
    // use dt from ATS for testing
    double dt = elm_coordinator_->get_dt(false);
    // call main method
    advance(&dt);
  }
}

void ELM_ATSDriver::finalize()
{
  WriteStateStatistics(*S_, *vo_);
  report_memory();
  Teuchos::TimeMonitor::summarize(*vo_->os());
  Coordinator::finalize();
}

} //namespace ATS

