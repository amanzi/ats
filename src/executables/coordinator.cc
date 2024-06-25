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
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "Units.hh"
#include "CompositeVector.hh"
#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "VisualizationDomainSet.hh"
#include "GeometricModel.hh"
#include "Checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PKFactory.hh"

#include "ats_mesh_factory.hh"

#include "coordinator.hh"

namespace ATS {

// this MUST be be called before using Coordinator
Coordinator::Coordinator(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                         const Amanzi::Comm_ptr_type& comm)
  : plist_(plist),
    comm_(comm),
    restart_(false),
    timer_(Teuchos::rcp(new Teuchos::Time("wallclock_monitor", true))),
    setup_timer_(Teuchos::TimeMonitor::getNewCounter("setup")),
    cycle_timer_(Teuchos::TimeMonitor::getNewCounter("cycle"))
{
  // create state.
  S_ = Teuchos::rcp(new Amanzi::State(Teuchos::sublist(plist_, "state")));

  // create the geometric model and regions
  Teuchos::ParameterList reg_list = plist_->sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_list, *comm_));

  // create and register meshes
  ATS::Mesh::createMeshes(plist_, comm_, gm, *S_);

  coordinator_list_ = Teuchos::sublist(plist_, "cycle driver");
  InitializeFromPlist_();

  // create the top level PK
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(plist_, "PKs");
  Teuchos::ParameterList pk_tree_list = coordinator_list_->sublist("PK tree");
  if (pk_tree_list.numParams() != 1) {
    Errors::Message message("CycleDriver: PK tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list.begin();
  const std::string& pk_name = pk_tree_list.name(pk_item);

  // create the pk
  Amanzi::PKFactory pk_factory;
  pk_ = pk_factory.CreatePK(pk_name, comm_, pk_tree_list, plist_, S_);
  pk_->modifyParameterList();
  pk_->parseParameterList();

  // create the checkpointing
  Teuchos::ParameterList& chkp_plist = plist_->sublist("checkpoint");
  checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(chkp_plist, *S_));

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

  for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
    if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
      Amanzi::Key node_key = Amanzi::Keys::getKey(mesh->first, "vertex_coordinates");
      S_->Require<Amanzi::CompositeVector, Amanzi::CompositeVectorSpace>(
          node_key, Amanzi::Tags::NEXT, node_key)
        .SetMesh(mesh->second.first)
        ->SetGhosted()
        ->SetComponent(
          "node", Amanzi::AmanziMesh::Entity_kind::NODE, mesh->second.first->getSpaceDimension());
    }

    // // writes region information
    // if (plist_->isSublist("analysis")) {
    //   Amanzi::InputAnalysis analysis(mesh->second.first, mesh->first);
    //   analysis.Init(plist_->sublist("analysis").sublist(mesh->first));
    //   analysis.RegionAnalysis();
    //   analysis.OutputBCs();
    // }
  }

  // create verbose object
  vo_ = Teuchos::rcp(new Amanzi::VerboseObject(comm_, "Coordinator", *coordinator_list_));

  // create the time step manager
  tsm_ = Teuchos::rcp(new Amanzi::TimeStepManager());
}

void
Coordinator::setup()
{
  // common constants
  S_->Require<double>("atmospheric_pressure", Amanzi::Tags::DEFAULT, "coordinator");
  S_->Require<Amanzi::AmanziGeometry::Point>("gravity", Amanzi::Tags::DEFAULT, "coordinator");

  // needed other times
  S_->require_time(Amanzi::Tags::CURRENT);
  S_->require_time(Amanzi::Tags::NEXT);

  // order matters here -- PKs set the leaves, then observations can use those
  // if provided, and setup finally deals with all secondaries and allocates memory
  pk_->setTags(Amanzi::Tags::CURRENT, Amanzi::Tags::NEXT);
  pk_->setup();
  for (auto& obs : observations_) obs->Setup(S_.ptr());
  S_->Setup();
}

void
Coordinator::initialize()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  int size = comm_->getSize();
  int rank = comm_->getRank();

  S_->set_time(Amanzi::Tags::CURRENT, t0_);
  S_->set_time(Amanzi::Tags::NEXT, t0_);
  S_->set_cycle(cycle0_);

  // Restart from checkpoint part 1:
  //  - get the time prior to initializing anything else
  if (restart_) {
    Amanzi::Checkpoint chkp(restart_filename_, comm_);
    double t_restart;
    chkp.read(Teuchos::ParameterList("time"), t_restart);
    S_->set_time(Amanzi::Tags::CURRENT, t_restart);
    S_->set_time(Amanzi::Tags::NEXT, t_restart);
    t0_ = t_restart;
  }

  // Initialize the state
  S_->InitializeFields();

  // Initialize the process kernels
  pk_->initialize();

  // calling CommitStep to set up copies as needed.
  pk_->commitStep(t0_, t0_, Amanzi::Tags::NEXT);

  // initialize vertex coordinate data
  for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
    if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
      Amanzi::Key node_key = Amanzi::Keys::getKey(mesh->first, "vertex_coordinates");
      copyMeshCoordinatesToVector(
        *mesh->second.first,
        Amanzi::AmanziMesh::Entity_kind::NODE,
        S_->GetW<Amanzi::CompositeVector>(node_key, Amanzi::Tags::NEXT, node_key));
      S_->GetRecordW(node_key, Amanzi::Tags::NEXT, node_key).set_initialized();
    }
  }

  // Restart from checkpoint part 2:
  // -- load all other data
  if (restart_) {
    Amanzi::Checkpoint chkp(restart_filename_, *S_);
    chkp.read(*S_);
    t0_ = S_->get_time();
    cycle0_ = S_->get_cycle();

    S_->set_time(Amanzi::Tags::CURRENT, t0_);
    S_->set_time(Amanzi::Tags::NEXT, t0_);

    for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
        Amanzi::Key node_key = Amanzi::Keys::getKey(mesh->first, "vertex_coordinates");
        copyVectorToMeshCoordinates(S_->Get<Amanzi::CompositeVector>(node_key, Amanzi::Tags::NEXT),
                                    *mesh->second.first);
      }
    }
  }

  // Final checks.
  //S_->CheckNotEvaluatedFieldsInitialized();
  S_->InitializeEvaluators();
  S_->InitializeFieldCopies();
  S_->CheckAllFieldsInitialized();

  // commit one more time, since some variables may have changed in the
  // previous call (test this... maybe chemistry/transport variables?  Is this
  // still necessary?  And do we need to set cycle to -1 here too? --ETC)
  pk_->commitStep(S_->get_time(), S_->get_time(), Amanzi::Tags::NEXT);

  // Write dependency graph.
  S_->WriteDependencyGraph();
  S_->InitializeIOFlags();

  // Check final initialization
  S_->WriteStatistics(vo_.ptr());

  // Set up visualization
  auto vis_list = Teuchos::sublist(plist_, "visualization");
  for (auto& entry : *vis_list) {
    std::string domain_name = entry.first;

    if (S_->HasMesh(domain_name)) {
      // visualize standard domain
      auto mesh_p = S_->GetMesh(domain_name);
      auto sublist_p = Teuchos::sublist(vis_list, domain_name);
      if (!sublist_p->isParameter("file name base")) {
        sublist_p->set<std::string>("file name base", "ats_vis");
      }

      if (S_->HasMesh(domain_name + "_3d") && sublist_p->get<bool>("visualize on 3D mesh", true))
        mesh_p = S_->GetMesh(domain_name + "_3d");

      // vis successful timesteps
      auto vis = Teuchos::rcp(new Amanzi::Visualization(*sublist_p, mesh_p, false));
      visualization_.push_back(vis);

    } else if (Amanzi::Keys::isDomainSet(domain_name)) {
      // visualize domain set
      const auto& dset = S_->GetDomainSet(Amanzi::Keys::getDomainSetName(domain_name));
      auto sublist_p = Teuchos::sublist(vis_list, domain_name);

      if (sublist_p->get("visualize individually", false)) {
        // visualize each subdomain
        for (const auto& subdomain : *dset) {
          Teuchos::ParameterList sublist = vis_list->sublist(subdomain);
          if (!sublist.isParameter("file name base"))
            sublist.set<std::string>("file name base", "ats_vis");
          auto vis =
            Teuchos::rcp(new Amanzi::Visualization(sublist, S_->GetMesh(subdomain), false));
          visualization_.push_back(vis);
        }
      } else {
        // visualize collectively
        auto domain_name_base = Amanzi::Keys::getDomainSetName(domain_name);
        if (!sublist_p->isParameter("file name base"))
          sublist_p->set<std::string>("file name base", "ats_vis");
        auto vis = Teuchos::rcp(new Amanzi::VisualizationDomainSet(
          *sublist_p, dset->getReferencingParent(), false, dset));
        visualization_.push_back(vis);
      }
    }
  }

  // make observations at time 0
  for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());

  // set up the TSM
  // -- register visualization times
  for (const auto& vis : visualization_) vis->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register checkpoint times
  checkpoint_->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register observation times
  for (const auto& obs : observations_) obs->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register the final time
  tsm_->RegisterTimeEvent(t1_);

  // -- register any intermediate requested times
  if (coordinator_list_->isSublist("required times")) {
    Teuchos::ParameterList& sublist = coordinator_list_->sublist("required times");
    Amanzi::IOEvent pause_times(sublist);
    pause_times.RegisterWithTimeStepManager(tsm_.ptr());
  }

  // -- advance cycle to 0 and begin
  if (S_->get_cycle() == -1) S_->advance_cycle();
}


void
Coordinator::finalize()
{
  // Force checkpoint at the end of simulation, and copy to checkpoint_final
  pk_->calculateDiagnostics(Amanzi::Tags::NEXT);
  checkpoint_->write(*S_, Amanzi::Checkpoint::WriteType::FINAL);

  // flush observations to make sure they are saved
  for (const auto& obs : observations_) obs->Flush();
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
      auto cell_map = (mesh->second.first)->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false);
      global_ncells += cell_map->getGlobalNumElements();
      local_ncells += cell_map->getLocalNumElements();
    }

    double mem = rss_usage();

    double percell(mem);
    if (local_ncells > 0) { percell = mem / local_ncells; }

    double max_percell(0.0);
    double min_percell(0.0);
    Teuchos::reduceAll(*comm_, Teuchos::REDUCE_MIN, 1, &percell, &min_percell);
    Teuchos::reduceAll(*comm_, Teuchos::REDUCE_MAX, 1, &percell, &max_percell);

    double total_mem(0.0);
    double max_mem(0.0);
    double min_mem(0.0);
    Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &mem, &total_mem);
    Teuchos::reduceAll(*comm_, Teuchos::REDUCE_MIN, 1, &mem, &min_mem);
    Teuchos::reduceAll(*comm_, Teuchos::REDUCE_MAX, 1, &mem, &max_mem);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "======================================================================"
               << std::endl;
    *vo_->os() << "All meshes combined have " << global_ncells << " cells." << std::endl;
    *vo_->os() << "Memory usage (high water mark):" << std::endl;
    *vo_->os() << std::fixed << std::setprecision(1);
    *vo_->os() << "  Maximum per core:   " << std::setw(7) << max_mem
               << " MBytes,  maximum per cell: " << std::setw(7) << max_percell * 1024 * 1024
               << " Bytes" << std::endl;
    *vo_->os() << "  Minimum per core:   " << std::setw(7) << min_mem
               << " MBytes,  minimum per cell: " << std::setw(7) << min_percell * 1024 * 1024
               << " Bytes" << std::endl;
    *vo_->os() << "  Total:              " << std::setw(7) << total_mem
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
Coordinator::InitializeFromPlist_()
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
Coordinator::getDt(bool after_fail)
{
  // get the physical step size
  double dt = pk_->getDt();
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
  double dt = S_->Get<double>("dt", Amanzi::Tags::DEFAULT);
  double t_old = S_->get_time(Amanzi::Tags::CURRENT);
  double t_new = S_->get_time(Amanzi::Tags::NEXT);

  bool fail = pk_->advanceStep(t_old, t_new, false);
  if (!fail) fail |= !pk_->isValidStep();

  // write state post-advance, if extreme
  S_->WriteStatistics(vo_.ptr(), Teuchos::VERB_EXTREME);

  if (!fail) {
    // commit the state, copying NEXT --> CURRENT
    pk_->commitStep(t_old, t_new, Amanzi::Tags::NEXT);

  } else {
    // Failed the timestep.
    // Potentially write out failed timestep for debugging
    for (const auto& vis : failed_visualization_) vis->write(*S_);

    // copy from old time into new time to reset the timestep
    pk_->failStep(t_old, t_new, Amanzi::Tags::NEXT);

    // check whether meshes are deformable, and if so, recover the old coordinates
    for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
        // collect the old coordinates
        Amanzi::Key node_key = Amanzi::Keys::getKey(mesh->first, "vertex_coordinates");
        Teuchos::RCP<const Amanzi::CompositeVector> vc_vec =
          S_->GetPtr<Amanzi::CompositeVector>(node_key, Amanzi::Tags::DEFAULT);
        vc_vec->scatterMasterToGhosted();
        copyVectorToMeshCoordinates(*vc_vec, *mesh->second.first);
      }
    }
  }
  // write state one extreme, post-commit/fail
  S_->WriteStatistics(vo_.ptr(), Teuchos::VERB_EXTREME);

  return fail;
}


bool
Coordinator::visualize(bool force)
{
  // write visualization if requested
  bool dump = force;
  int cycle = S_->get_cycle();
  double time = S_->get_time();

  if (!dump) {
    for (const auto& vis : visualization_) { dump |= vis->DumpRequested(cycle, time); }
  }

  if (dump) { pk_->calculateDiagnostics(Amanzi::Tags::NEXT); }

  for (const auto& vis : visualization_) {
    if (force || vis->DumpRequested(cycle, time)) { vis->write(*S_); }
  }
  return dump;
}


bool
Coordinator::checkpoint(bool force)
{
  int cycle = S_->get_cycle();
  double time = S_->get_time();
  bool dump = force;
  dump |= checkpoint_->DumpRequested(cycle, time);
  if (dump) checkpoint_->write(*S_);
  return dump;
}

} // namespace ATS
