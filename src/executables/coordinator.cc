/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Implementation for the Coordinator. Coordinator holds the functionality
called by the cycle driver, which runs the overall, top level timestep loop.

Coordinator instantiates states, ensures they are initialized, advances
timesteps, and writes vis and restart/checkpoint dumps. It contains one and
only one PK -- most likely this PK is an MPC of some type -- to do the
actual work.
------------------------------------------------------------------------- */

#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include "errors.hh"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Teuchos_CommHelpers.hpp"

#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "InputAnalysis.hh"

#include "Units.hh"
#include "CompositeVector.hh"
#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "VisualizationDomainSet.hh"
#include "IO.hh"
#include "GeometricModel.hh"
#include "Checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"
#include "pk_helpers.hh"

#include "ats_mesh_factory.hh"

#include "coordinator.hh"

namespace ATS {

// this MUST be be called before using Coordinator
Coordinator::Coordinator(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                         const Teuchos::RCP<Teuchos::Time>& wallclock_timer,
                         const Teuchos::RCP<const Teuchos::Comm<int>>& teuchos_comm,
                         const Amanzi::Comm_ptr_type& comm)
  : plist_(plist),
    teuchos_comm_(teuchos_comm),
    comm_(comm),
    restart_(false),
    wallclock_monitor_(*wallclock_timer),
    wallclock_timer_(wallclock_timer)
{
  // create verbose object and timers
  coordinator_list_ = Teuchos::sublist(plist_, "cycle driver");
  vo_ = Teuchos::rcp(new Amanzi::VerboseObject(comm_, "ATS", *coordinator_list_));
  Teuchos::OSTab tab = vo_->getOSTab();

  timers_["0: create mesh"] = Teuchos::TimeMonitor::getNewCounter("0: create mesh");
  timers_["1: create run"] = Teuchos::TimeMonitor::getNewCounter("1: create run");
  timers_["2: setup"] = Teuchos::TimeMonitor::getNewCounter("2: setup");
  timers_["3: initialize"] = Teuchos::TimeMonitor::getNewCounter("3: initialize");
  timers_["4: solve"] = Teuchos::TimeMonitor::getNewCounter("4: solve");
  timers_["4a: advance step"] = Teuchos::TimeMonitor::getNewCounter("4a: advance step");
  timers_["4b: visualize"] = Teuchos::TimeMonitor::getNewCounter("4b: visualize");
  timers_["4c: observe"] = Teuchos::TimeMonitor::getNewCounter("4c: observe");
  timers_["4d: checkpoint"] = Teuchos::TimeMonitor::getNewCounter("4d: checkpoint");
  timers_["5: finalize"] = Teuchos::TimeMonitor::getNewCounter("5: finalize");

  // print header material
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "Writing input file ..." << std::endl << std::endl;
    Teuchos::writeParameterListToXmlOStream(*plist_, *vo_->os());
    *vo_->os() << "  ... completed." << std::endl;
  }

  // construct state, geometric model, meshes
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning mesh creation stage..." << std::endl
               << std::flush;
  }
  {
    Teuchos::TimeMonitor timer(*timers_.at("0: create mesh"));

    // create state.
    S_ = Teuchos::rcp(new Amanzi::State(plist_->sublist("state")));

    // create the geometric model and regions
    Teuchos::ParameterList reg_list = plist_->sublist("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_list, *comm_));

    // create and register meshes
    ATS::Mesh::createMeshes(*plist_, comm_, gm, *S_);
  }
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("0: create mesh");
  }


  // create PKs, etc
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning run creation stage..." << std::endl
               << std::flush;
  }
  {
    Teuchos::TimeMonitor timer(*timers_.at("1: create run"));
    initializeFromPlist_();

    // create the top level PK
    Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(plist_, "PKs");
    Teuchos::ParameterList pk_tree_list = coordinator_list_->sublist("PK tree");
    if (pk_tree_list.numParams() != 1) {
      Errors::Message message(
        "CycleDriver: PK tree list should contain exactly one root node list");
      Exceptions::amanzi_throw(message);
    }
    Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list.begin();
    const std::string& pk_name = pk_tree_list.name(pk_item);

    // create the solution
    soln_ = Teuchos::rcp(new Amanzi::TreeVector(comm_));

    // create the pk
    Amanzi::PKFactory pk_factory;
    pk_ = pk_factory.CreatePK(pk_name, pk_tree_list, plist_, S_, soln_);

    // create the checkpointing
    Teuchos::ParameterList& chkp_plist = plist_->sublist("checkpoint");
    checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(chkp_plist, *S_));

    for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
        Amanzi::Key node_key = Amanzi::Keys::getKey(mesh->first, "vertex_coordinates");
        S_->Require<Amanzi::CompositeVector, Amanzi::CompositeVectorSpace>(
            node_key, Amanzi::Tags::NEXT, node_key)
          .SetMesh(mesh->second.first)
          ->SetGhosted()
          ->SetComponent("node", Amanzi::AmanziMesh::NODE, mesh->second.first->space_dimension());
      }

      // writes region information
      if (plist_->isSublist("analysis")) {
        Amanzi::InputAnalysis analysis(mesh->second.first, mesh->first);
        analysis.Init(plist_->sublist("analysis").sublist(mesh->first));
        analysis.RegionAnalysis();
        analysis.OutputBCs();
      }
    }

    // create the observations
    Teuchos::ParameterList& observation_plist = plist_->sublist("observations");
    for (auto& sublist : observation_plist) {
      if (observation_plist.isSublist(sublist.first)) {
        observations_.emplace_back(Teuchos::rcp(
          new Amanzi::UnstructuredObservations(observation_plist.sublist(sublist.first))));
      } else {
        Errors::Message msg("\"observations\" list must only include sublists.");
        Exceptions::amanzi_throw(msg);
      }
    }

    // create visualization
    auto vis_list = Teuchos::sublist(plist_, "visualization");
    for (auto& entry : *vis_list) {
      std::string domain_name = entry.first;

      if (S_->HasMesh(domain_name)) {
        // visualize standard domain
        auto mesh_p = S_->GetMesh(domain_name);
        auto sublist_p = Teuchos::sublist(vis_list, domain_name);
        if (!sublist_p->isParameter("file name base")) {
          if (domain_name.empty() || domain_name == "domain") {
            sublist_p->set<std::string>("file name base", std::string("ats_vis"));
          } else {
            sublist_p->set<std::string>("file name base", std::string("ats_vis_") + domain_name);
          }
        }

        if (S_->HasMesh(domain_name + "_3d") && sublist_p->get<bool>("visualize on 3D mesh", true))
          mesh_p = S_->GetMesh(domain_name + "_3d");

        // vis successful timesteps
        auto vis = Teuchos::rcp(new Amanzi::Visualization(*sublist_p));
        vis->set_name(domain_name);
        vis->set_mesh(mesh_p);
        vis->CreateFiles(false);
        visualization_.push_back(vis);

      } else if (Amanzi::Keys::isDomainSet(domain_name)) {
        // visualize domain set
        const auto& dset = S_->GetDomainSet(Amanzi::Keys::getDomainSetName(domain_name));
        auto sublist_p = Teuchos::sublist(vis_list, domain_name);

        if (sublist_p->get("visualize individually", false)) {
          // visualize each subdomain
          for (const auto& subdomain : *dset) {
            Teuchos::ParameterList sublist = vis_list->sublist(subdomain);
            sublist.set<std::string>("file name base", std::string("ats_vis_") + subdomain);
            auto vis = Teuchos::rcp(new Amanzi::Visualization(sublist));
            vis->set_name(subdomain);
            vis->set_mesh(S_->GetMesh(subdomain));
            vis->CreateFiles(false);
            visualization_.push_back(vis);
          }
        } else {
          // visualize collectively
          auto domain_name_base = Amanzi::Keys::getDomainSetName(domain_name);
          if (!sublist_p->isParameter("file name base"))
            sublist_p->set("file name base", std::string("ats_vis_") + domain_name_base);
          auto vis = Teuchos::rcp(new Amanzi::VisualizationDomainSet(*sublist_p));
          vis->set_name(domain_name_base);
          vis->set_domain_set(dset);
          vis->set_mesh(dset->get_referencing_parent());
          vis->CreateFiles(false);
          visualization_.push_back(vis);
        }
      }
    }

    // create the time step manager, register assorted timestep control events
    tsm_ = Teuchos::rcp(new Amanzi::TimeStepManager());
    checkpoint_->RegisterWithTimeStepManager(tsm_.ptr());
    for (const auto& obs : observations_) obs->RegisterWithTimeStepManager(tsm_.ptr());
    for (const auto& vis : visualization_) vis->RegisterWithTimeStepManager(tsm_.ptr());
    tsm_->RegisterTimeEvent(t1_);

    if (coordinator_list_->isSublist("required times")) {
      Teuchos::ParameterList& sublist = coordinator_list_->sublist("required times");
      Amanzi::IOEvent pause_times(sublist);
      pause_times.RegisterWithTimeStepManager(tsm_.ptr());
    }
  }
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("1: create run");
  }
}


void
Coordinator::setup()
{
  Teuchos::TimeMonitor monitor(*timers_.at("2: setup"));

  // common constants
  S_->Require<double>("atmospheric_pressure", Amanzi::Tags::DEFAULT, "coordinator");
  S_->Require<Amanzi::AmanziGeometry::Point>("gravity", Amanzi::Tags::DEFAULT, "coordinator");

  // needed other times
  S_->require_time(Amanzi::Tags::CURRENT);
  S_->require_time(Amanzi::Tags::NEXT);

  // order matters here -- PKs set the leaves, then observations can use those
  // if provided, and setup finally deals with all secondaries and allocates memory
  pk_->set_tags(Amanzi::Tags::CURRENT, Amanzi::Tags::NEXT);
  pk_->Setup();
  for (auto& obs : observations_) obs->Setup(S_.ptr());
  S_->Setup();
}

void
Coordinator::initialize()
{
  Teuchos::TimeMonitor monitor(*timers_.at("3: initialize"));
  Teuchos::OSTab tab = vo_->getOSTab();

  int size = comm_->NumProc();
  int rank = comm_->MyPID();

  S_->set_time(Amanzi::Tags::CURRENT, t0_);
  S_->set_time(Amanzi::Tags::NEXT, t0_);
  S_->set_cycle(cycle0_);

  // Restart from checkpoint part 1:
  //  - get the time prior to initializing anything else
  if (restart_) {
    double t_restart = Amanzi::ReadCheckpointInitialTime(comm_, restart_filename_);
    S_->set_time(Amanzi::Tags::CURRENT, t_restart);
    S_->set_time(Amanzi::Tags::NEXT, t_restart);
    t0_ = t_restart;
  }

  // Initialize the state
  S_->InitializeFields();

  // Initialize the process kernels
  pk_->Initialize();

  // calling CommitStep to set up copies as needed.
  pk_->CommitStep(t0_, t0_, Amanzi::Tags::NEXT);

  // initialize vertex coordinate data
  for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
    if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
      Amanzi::Key node_key = Amanzi::Keys::getKey(mesh->first, "vertex_coordinates");
      copyMeshCoordinatesToVector(
        *mesh->second.first,
        S_->GetW<Amanzi::CompositeVector>(node_key, Amanzi::Tags::NEXT, node_key));
      S_->GetRecordW(node_key, Amanzi::Tags::NEXT, node_key).set_initialized();
    }
  }

  // Restart from checkpoint part 2:
  // -- load all other data
  if (restart_) {
    Amanzi::ReadCheckpoint(comm_, *S_, restart_filename_);
    t0_ = S_->get_time();
    cycle0_ = S_->get_cycle();

    S_->set_time(Amanzi::Tags::CURRENT, t0_);
    S_->set_time(Amanzi::Tags::NEXT, t0_);

    for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first))
        Amanzi::DeformCheckpointMesh(*S_, mesh->first);
    }
  }

  // Final checks.
  //S_->CheckNotEvaluatedFieldsInitialized();
  S_->InitializeEvaluators();
  S_->InitializeFieldCopies();
  S_->CheckAllFieldsInitialized();
  S_->InitializeIOFlags();

  // commit one more time, since some variables may have changed in the
  // previous call (test this... maybe chemistry/transport variables?  Is this
  // still necessary?  And do we need to set cycle to -1 here too? --ETC)
  pk_->CommitStep(S_->get_time(), S_->get_time(), Amanzi::Tags::NEXT);

  // -- advance cycle to 0 and begin
  if (S_->get_cycle() == -1) S_->advance_cycle();
}


void
Coordinator::finalize()
{
  Teuchos::TimeMonitor monitor(*timers_.at("5: finalize"));

  // Force checkpoint at the end of simulation, and copy to checkpoint_final
  pk_->CalculateDiagnostics(Amanzi::Tags::NEXT);
  checkpoint_->Write(*S_, Amanzi::Checkpoint::WriteType::FINAL);

  // flush observations to make sure they are saved
  for (const auto& obs : observations_) obs->Flush();

  // report out
  WriteStateStatistics(*S_, *vo_);
  report_memory();
}


double
rss_usage()
{ // return ru_maxrss in MBytes
#if (defined(__unix__) || defined(__unix) || defined(unix) || defined(__APPLE__) ||                \
     defined(__MACH__))
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#  if (defined(__APPLE__) || defined(__MACH__))
  return static_cast<double>(usage.ru_maxrss) / 1024.0 / 1024.0;
#  else
  return static_cast<double>(usage.ru_maxrss) / 1024.0;
#  endif
#else
  return 0.0;
#endif
}


void
Coordinator::report_memory()
{
  // report the memory high water mark (using ru_maxrss)
  // this should be called at the very end of a simulation
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    double global_ncells(0.0);
    double local_ncells(0.0);
    for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
      Epetra_Map cell_map = (mesh->second.first)->cell_map(false);
      global_ncells += cell_map.NumGlobalElements();
      local_ncells += cell_map.NumMyElements();
    }

    double mem = rss_usage();

    double percell(mem);
    if (local_ncells > 0) { percell = mem / local_ncells; }

    double max_percell(0.0);
    double min_percell(0.0);
    comm_->MinAll(&percell, &min_percell, 1);
    comm_->MaxAll(&percell, &max_percell, 1);

    double total_mem(0.0);
    double max_mem(0.0);
    double min_mem(0.0);
    comm_->SumAll(&mem, &total_mem, 1);
    comm_->MinAll(&mem, &min_mem, 1);
    comm_->MaxAll(&mem, &max_mem, 1);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "--------------------------------------------------------------------------------"
               << std::endl
               << "All meshes combined have " << global_ncells << " cells." << std::endl
               << "Memory usage (high water mark):" << std::endl
               << std::fixed << std::setprecision(1) << "  Maximum per core:   " << std::setw(7)
               << max_mem << " MBytes,  maximum per cell: " << std::setw(7)
               << max_percell * 1024 * 1024 << " Bytes" << std::endl
               << "  Minimum per core:   " << std::setw(7) << min_mem
               << " MBytes,  minimum per cell: " << std::setw(7) << min_percell * 1024 * 1024
               << " Bytes" << std::endl
               << "  Total:              " << std::setw(7) << total_mem
               << " MBytes,  total per cell:   " << std::setw(7)
               << total_mem / global_ncells * 1024 * 1024 << " Bytes" << std::endl;
  }


  // double doubles_count(0.0);
  // for (Amanzi::State::data_iterator field=S_->data_begin(); field!=S_->data_end(); ++field) {
  //   doubles_count += static_cast<double>(field->second->GetLocalElementCount());
  // }
  // double global_doubles_count(0.0);
  // double min_doubles_count(0.0);
  // double max_doubles_count(0.0);
  // comm_->SumAll(&doubles_count,&global_doubles_count,1);
  // comm_->MinAll(&doubles_count,&min_doubles_count,1);
  // comm_->MaxAll(&doubles_count,&max_doubles_count,1);

  // Teuchos::OSTab tab = vo_->getOSTab();
  // *vo_->os() << "Doubles allocated in state fields " << std::endl;
  // *vo_->os() << "  Maximum per core:   " << std::setw(7)
  //            << max_doubles_count*8/1024/1024 << " MBytes" << std::endl;
  // *vo_->os() << "  Minimum per core:   " << std::setw(7)
  //            << min_doubles_count*8/1024/1024 << " MBytes" << std::endl;
  // *vo_->os() << "  Total:              " << std::setw(7)
  //            << global_doubles_count*8/1024/1024 << " MBytes" << std::endl;
}


void
Coordinator::initializeFromPlist_()
{
  Amanzi::Utils::Units units;
  t0_ = coordinator_list_->get<double>("start time");
  std::string t0_units = coordinator_list_->get<std::string>("start time units", "s");
  if (!units.IsValidTime(t0_units)) {
    Errors::Message msg;
    msg << "Coordinator start time: unknown time units type: \"" << t0_units
        << "\"  Valid are: " << units.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
  bool success;
  t0_ = units.ConvertTime(t0_, t0_units, "s", success);

  t1_ = coordinator_list_->get<double>("end time");
  std::string t1_units = coordinator_list_->get<std::string>("end time units", "s");
  if (!units.IsValidTime(t1_units)) {
    Errors::Message msg;
    msg << "Coordinator end time: unknown time units type: \"" << t1_units
        << "\"  Valid are: " << units.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
  t1_ = units.ConvertTime(t1_, t1_units, "s", success);

  max_dt_ = coordinator_list_->get<double>("max time step size [s]", 1.0e99);
  min_dt_ = coordinator_list_->get<double>("min time step size [s]", 1.0e-12);
  cycle0_ = coordinator_list_->get<int>("start cycle", -1);
  cycle1_ = coordinator_list_->get<int>("end cycle", -1);
  duration_ = coordinator_list_->get<double>("wallclock duration [hrs]", -1.0);
  subcycled_ts_ = coordinator_list_->get<bool>("subcycled timestep", false);

  // restart control
  restart_ = coordinator_list_->isParameter("restart from checkpoint file");
  if (restart_)
    restart_filename_ = coordinator_list_->get<std::string>("restart from checkpoint file");
}


// -----------------------------------------------------------------------------
// acquire the chosen timestep size
// -----------------------------------------------------------------------------
double
Coordinator::get_dt(bool after_fail)
{
  // get the physical step size
  double dt = pk_->get_dt();
  double dt_pk = dt;
  if (dt < 0.) return dt;

  // check if the step size has gotten too small
  if (dt < min_dt_) {
    Errors::Message message("Coordinator: error, timestep too small");
    Exceptions::amanzi_throw(message);
  }

  // cap the max step size
  if (dt > max_dt_) { dt = max_dt_; }

  // ask the step manager if this step is ok
  dt = tsm_->TimeStep(S_->get_time(Amanzi::Tags::NEXT), dt, after_fail);

  // note, I believe this can go away (along with the input spec flag) once
  // amanzi/amanzi#685 is closed --etc
  if (subcycled_ts_) dt = std::min(dt, dt_pk);
  return dt;
}


bool
Coordinator::advance()
{
  Teuchos::TimeMonitor timer(*timers_.at("4a: advance step"));
  double dt = S_->Get<double>("dt", Amanzi::Tags::DEFAULT);
  double t_old = S_->get_time(Amanzi::Tags::CURRENT);
  double t_new = S_->get_time(Amanzi::Tags::NEXT);

  bool fail = pk_->AdvanceStep(t_old, t_new, false);
  if (!fail) fail |= !pk_->ValidStep();

  // write state post-advance, if extreme
  WriteStateStatistics(*S_, *vo_, Teuchos::VERB_EXTREME);

  if (!fail) {
    // commit the state, copying NEXT --> CURRENT
    pk_->CommitStep(t_old, t_new, Amanzi::Tags::NEXT);

  } else {
    // Failed the timestep.
    // Potentially write out failed timestep for debugging
    for (const auto& vis : failed_visualization_) WriteVis(*vis, *S_);

    // copy from old time into new time to reset the timestep
    pk_->FailStep(t_old, t_new, Amanzi::Tags::NEXT);

    // check whether meshes are deformable, and if so, recover the old coordinates
    for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
        // collect the old coordinates
        Amanzi::Key node_key = Amanzi::Keys::getKey(mesh->first, "vertex_coordinates");
        Teuchos::RCP<const Amanzi::CompositeVector> vc_vec =
          S_->GetPtr<Amanzi::CompositeVector>(node_key, Amanzi::Tags::DEFAULT);
        vc_vec->ScatterMasterToGhosted();
        const Epetra_MultiVector& vc = *vc_vec->ViewComponent("node", true);

        std::vector<int> node_ids(vc.MyLength());
        Amanzi::AmanziGeometry::Point_List old_positions(vc.MyLength());
        for (int n = 0; n != vc.MyLength(); ++n) {
          node_ids[n] = n;
          if (mesh->second.first->space_dimension() == 2) {
            old_positions[n] = Amanzi::AmanziGeometry::Point(vc[0][n], vc[1][n]);
          } else {
            old_positions[n] = Amanzi::AmanziGeometry::Point(vc[0][n], vc[1][n], vc[2][n]);
          }
        }

        // undeform the mesh
        Amanzi::AmanziGeometry::Point_List final_positions;
        mesh->second.first->deform(node_ids, old_positions, false, &final_positions);
      }
    }
  }
  // write state one extreme, post-commit/fail
  WriteStateStatistics(*S_, *vo_, Teuchos::VERB_EXTREME);

  return fail;
}


bool
Coordinator::visualize(bool force)
{
  Teuchos::TimeMonitor timer(*timers_.at("4b: visualize"));
  // write visualization if requested
  bool dump = force;
  int cycle = S_->get_cycle();
  double time = S_->get_time();

  if (!dump) {
    for (const auto& vis : visualization_) { dump |= vis->DumpRequested(cycle, time); }
  }

  if (dump) { pk_->CalculateDiagnostics(Amanzi::Tags::NEXT); }

  for (const auto& vis : visualization_) {
    if (force || vis->DumpRequested(cycle, time)) { WriteVis(*vis, *S_); }
  }
  return dump;
}


void
Coordinator::observe()
{
  Teuchos::TimeMonitor timer(*timers_.at("4c: observe"));
  // make observations, vis, and checkpoints
  for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());
}

bool
Coordinator::checkpoint(bool force)
{
  Teuchos::TimeMonitor timer(*timers_.at("4d: checkpoint"));
  int cycle = S_->get_cycle();
  double time = S_->get_time();
  bool dump = force;
  dump |= checkpoint_->DumpRequested(cycle, time);
  if (dump) checkpoint_->Write(*S_);
  return dump;
}

void
Coordinator::reportOneTimer_(const std::string& timer_name)
{
  auto timer = timers_.at(timer_name);
  double l_time = timer->totalElapsedTime();
  double min_time, max_time, mean_time;
  Teuchos::reduceAll(*teuchos_comm_, Teuchos::REDUCE_MIN, 1, &l_time, &min_time);
  Teuchos::reduceAll(*teuchos_comm_, Teuchos::REDUCE_MAX, 1, &l_time, &max_time);
  Teuchos::reduceAll(*teuchos_comm_, Teuchos::REDUCE_SUM, 1, &l_time, &mean_time);
  mean_time = mean_time / teuchos_comm_->getSize();
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << min_time << " / " << mean_time << " / " << max_time << " (min/mean/max) [s]"
               << std::endl;
  }
}

} // namespace ATS
