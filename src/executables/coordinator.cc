/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the Coordinator.  Coordinator is basically just a class to hold
the cycle driver, which runs the overall, top level timestep loop.  It
instantiates states, ensures they are initialized, and runs the timestep loop
including Vis and restart/checkpoint dumps.  It contains one and only one PK
-- most likely this PK is an MPC of some type -- to do the actual work.
------------------------------------------------------------------------- */

#include <iostream>
#include <unistd.h>
#include <sys/resource.h>
#include "errors.hh"

#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "InputAnalysis.hh"

#include "Units.hh"

#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "VisualizationDomainSet.hh"
#include "Checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"

#include "coordinator.hh"

#define DEBUG_MODE 1

namespace ATS {

Coordinator::Coordinator(Teuchos::ParameterList& parameter_list,
                         Teuchos::RCP<Amanzi::State>& S,
                         Amanzi::Comm_ptr_type comm ) :
    parameter_list_(Teuchos::rcp(new Teuchos::ParameterList(parameter_list))),
    S_(S),
    comm_(comm),
    restart_(false) {

  // create and start the global timer
  timer_ = Teuchos::rcp(new Teuchos::Time("wallclock_monitor",true));
  setup_timer_ = Teuchos::TimeMonitor::getNewCounter("setup");
  cycle_timer_ = Teuchos::TimeMonitor::getNewCounter("cycle");
  coordinator_init();

  vo_ = Teuchos::rcp(new Amanzi::VerboseObject("Coordinator", *parameter_list_));
};

void Coordinator::coordinator_init() {
  coordinator_list_ = Teuchos::sublist(parameter_list_, "cycle driver");
  read_parameter_list();

  // create the top level PK
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(parameter_list_, "PKs");
  Teuchos::ParameterList pk_tree_list = coordinator_list_->sublist("PK tree");
  if (pk_tree_list.numParams() != 1) {
    Errors::Message message("CycleDriver: PK tree list should contain exactly one root node list");
    Exceptions::amanzi_throw(message);
  }
  Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list.begin();
  const std::string &pk_name = pk_tree_list.name(pk_item);

  // create the solution
  soln_ = Teuchos::rcp(new Amanzi::TreeVector());

  // create the pk
  Amanzi::PKFactory pk_factory;
  pk_ = pk_factory.CreatePK(pk_name, pk_tree_list, parameter_list_, S_, soln_);

  // create the checkpointing
  Teuchos::ParameterList& chkp_plist = parameter_list_->sublist("checkpoint");
  checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(chkp_plist, *S_));

  // create the observations
  Teuchos::ParameterList& observation_plist = parameter_list_->sublist("observations");
  for (auto& sublist : observation_plist) {
    if (observation_plist.isSublist(sublist.first)) {
      observations_.emplace_back(Teuchos::rcp(new Amanzi::UnstructuredObservations(
                observation_plist.sublist(sublist.first))));
    } else {
      Errors::Message msg("\"observations\" list must only include sublists.");
      Exceptions::amanzi_throw(msg);
    }
  }

  // check whether meshes are deformable, and if so require a nodal position
  for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
       mesh!=S_->mesh_end(); ++mesh) {

    if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
      std::string node_key;
      if (mesh->first != "domain") node_key= mesh->first+std::string("-vertex_coordinate");
      else node_key = std::string("vertex_coordinate");

      S_->RequireField(node_key)->SetMesh(mesh->second.first)->SetGhosted()
          ->AddComponent("node", Amanzi::AmanziMesh::NODE, mesh->second.first->space_dimension());
    }

    // writes region information
    if (parameter_list_->isSublist("analysis")) {
      Amanzi::InputAnalysis analysis(mesh->second.first, mesh->first);
      analysis.Init(parameter_list_->sublist("analysis").sublist(mesh->first));
      analysis.RegionAnalysis();
      analysis.OutputBCs();
    }
  }

  // create the time step manager
  tsm_ = Teuchos::rcp(new Amanzi::TimeStepManager());
}

void Coordinator::setup() {
  // Set up the states, creating all data structures.
  S_->set_time(t0_);
  S_->set_cycle(cycle0_);
  S_->RequireScalar("dt", "coordinator");

  // order matters here -- PKs set the leaves, then observations can use those
  // if provided, and setup finally deals with all secondaries and allocates memory
  pk_->Setup(S_.ptr());
  for (auto& obs : observations_) obs->Setup(S_.ptr());
  S_->Setup();
}

double Coordinator::initialize() {
  Teuchos::OSTab tab = vo_->getOSTab();
  int size = comm_->NumProc();
  int rank = comm_->MyPID();

  // Restart from checkpoint part 1:
  //  - get the time prior to initializing anything else
  if (restart_) {
    S_->set_time(Amanzi::ReadCheckpointInitialTime(comm_, restart_filename_));
  }

  // Initialize the state
  *S_->GetScalarData("dt", "coordinator") = 0.;
  S_->GetField("dt","coordinator")->set_initialized();
  S_->InitializeFields();

  // Initialize the process kernels
  pk_->Initialize(S_.ptr());

  // Restart from checkpoint part 2:
  // -- load all other data
  double dt_restart = -1;
  if (restart_) {
    dt_restart = Amanzi::ReadCheckpoint(*S_, restart_filename_);
    t0_ = S_->time();
    cycle0_ = S_->cycle();

    for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
         mesh!=S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first)) {
        Amanzi::DeformCheckpointMesh(*S_, mesh->first);
      }
    }
  }

  // Final checks.
  S_->CheckNotEvaluatedFieldsInitialized();
  S_->InitializeEvaluators();
  S_->InitializeFieldCopies();
  S_->CheckAllFieldsInitialized();

  // commit the initial conditions.
  pk_->CommitStep(S_->time(), S_->time(), S_);

  // Write dependency graph.
  S_->WriteDependencyGraph();
  S_->InitializeIOFlags();

  // Check final initialization
  WriteStateStatistics(*S_, *vo_);

  // Set up visualization
  auto vis_list = Teuchos::sublist(parameter_list_,"visualization");
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
          sublist_p->set<std::string>("file name base", std::string("ats_vis_")+domain_name);
        }
      }

      if (S_->HasMesh(domain_name+"_3d") && sublist_p->get<bool>("visualize on 3D mesh", true))
        mesh_p = S_->GetMesh(domain_name+"_3d");

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
          sublist.set<std::string>("file name base", std::string("ats_vis_")+subdomain);
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
          sublist_p->set("file name base", std::string("ats_vis_")+domain_name_base);
        auto vis = Teuchos::rcp(new Amanzi::VisualizationDomainSet(*sublist_p));
        vis->set_name(domain_name_base);
        vis->set_mesh(dset->get_referencing_parent());
        for (const auto& subdomain : *dset) {
          vis->set_subdomain_mesh(subdomain, S_->GetMesh(subdomain));
        }
        vis->CreateFiles(false);
        visualization_.push_back(vis);
      }
    }
  }

  // make observations at time 0
  for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());

  S_->set_time(t0_); // in case steady state solve changed this
  S_->set_cycle(cycle0_);

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

  // Create an intermediate state that will store the updated solution until
  // we know it has succeeded.
  S_next_ = Teuchos::rcp(new Amanzi::State(*S_));
  *S_next_ = *S_;
  if (subcycled_ts_) {
    S_inter_ = Teuchos::rcp(new Amanzi::State(*S_));
    *S_inter_ = *S_;
  } else {
    S_inter_ = S_;
  }

  // set the states in the PKs Passing null for S_ allows for safer subcycling
  // -- PKs can't use it, so it is guaranteed to be pristinely the old
  // timestep.  This code is useful for testing that PKs don't use S.
  // pk_->set_states(Teuchos::null, S_inter_, S_next_);

  // That said, S is required to be valid, since it is valid for all time
  // (e.g. from construction), whereas S_inter and S_next are only valid after
  // set_states() is called.  This allows for standard interfaces.
  pk_->set_states(S_, S_inter_, S_next_);
  return dt_restart;
}

void Coordinator::finalize() {
  // Force checkpoint at the end of simulation, and copy to checkpoint_final
  pk_->CalculateDiagnostics(S_next_);
  checkpoint_->Write(*S_next_, 0.0, true);

  // flush observations to make sure they are saved
  for (const auto& obs : observations_) obs->Flush();
}


double rss_usage() { // return ru_maxrss in MBytes
#if (defined(__unix__) || defined(__unix) || defined(unix) || defined(__APPLE__) || defined(__MACH__))
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#if (defined(__APPLE__) || defined(__MACH__))
  return static_cast<double>(usage.ru_maxrss)/1024.0/1024.0;
#else
  return static_cast<double>(usage.ru_maxrss)/1024.0;
#endif
#else
  return 0.0;
#endif
}


void Coordinator::report_memory() {
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
    if (local_ncells > 0) {
      percell = mem/local_ncells;
    }

    double max_percell(0.0);
    double min_percell(0.0);
    comm_->MinAll(&percell,&min_percell,1);
    comm_->MaxAll(&percell,&max_percell,1);

    double total_mem(0.0);
    double max_mem(0.0);
    double min_mem(0.0);
    comm_->SumAll(&mem,&total_mem,1);
    comm_->MinAll(&mem,&min_mem,1);
    comm_->MaxAll(&mem,&max_mem,1);

    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "======================================================================" << std::endl;
    *vo_->os() << "All meshes combined have " << global_ncells << " cells." << std::endl;
    *vo_->os() << "Memory usage (high water mark):" << std::endl;
    *vo_->os() << std::fixed << std::setprecision(1);
    *vo_->os() << "  Maximum per core:   " << std::setw(7) << max_mem
          << " MBytes,  maximum per cell: " << std::setw(7) << max_percell*1024*1024
          << " Bytes" << std::endl;
    *vo_->os() << "  Minimum per core:   " << std::setw(7) << min_mem
          << " MBytes,  minimum per cell: " << std::setw(7) << min_percell*1024*1024
         << " Bytes" << std::endl;
    *vo_->os() << "  Total:              " << std::setw(7) << total_mem
          << " MBytes,  total per cell:   " << std::setw(7) << total_mem/global_ncells*1024*1024
          << " Bytes" << std::endl;
  }


  double doubles_count(0.0);
  for (Amanzi::State::field_iterator field=S_->field_begin(); field!=S_->field_end(); ++field) {
    doubles_count += static_cast<double>(field->second->GetLocalElementCount());
  }
  double global_doubles_count(0.0);
  double min_doubles_count(0.0);
  double max_doubles_count(0.0);
  comm_->SumAll(&doubles_count,&global_doubles_count,1);
  comm_->MinAll(&doubles_count,&min_doubles_count,1);
  comm_->MaxAll(&doubles_count,&max_doubles_count,1);

  Teuchos::OSTab tab = vo_->getOSTab();
  *vo_->os() << "Doubles allocated in state fields " << std::endl;
  *vo_->os() << "  Maximum per core:   " << std::setw(7)
             << max_doubles_count*8/1024/1024 << " MBytes" << std::endl;
  *vo_->os() << "  Minimum per core:   " << std::setw(7)
             << min_doubles_count*8/1024/1024 << " MBytes" << std::endl;
  *vo_->os() << "  Total:              " << std::setw(7)
             << global_doubles_count*8/1024/1024 << " MBytes" << std::endl;
}



void Coordinator::read_parameter_list() {
  Amanzi::Utils::Units units;
  t0_ = coordinator_list_->get<double>("start time");
  std::string t0_units = coordinator_list_->get<std::string>("start time units", "s");
  if (!units.IsValidTime(t0_units)) {
    Errors::Message msg;
    msg << "Coordinator start time: unknown time units type: \"" << t0_units << "\"  Valid are: " << units.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
  bool success;
  t0_ = units.ConvertTime(t0_, t0_units, "s", success);

  t1_ = coordinator_list_->get<double>("end time");
  std::string t1_units = coordinator_list_->get<std::string>("end time units", "s");
  if (!units.IsValidTime(t1_units)) {
    Errors::Message msg;
    msg << "Coordinator end time: unknown time units type: \"" << t1_units << "\"  Valid are: " << units.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
  t1_ = units.ConvertTime(t1_, t1_units, "s", success);

  max_dt_ = coordinator_list_->get<double>("max time step size [s]", 1.0e99);
  min_dt_ = coordinator_list_->get<double>("min time step size [s]", 1.0e-12);
  cycle0_ = coordinator_list_->get<int>("start cycle",0);
  cycle1_ = coordinator_list_->get<int>("end cycle",-1);
  duration_ = coordinator_list_->get<double>("wallclock duration [hrs]", -1.0);
  subcycled_ts_ = coordinator_list_->get<bool>("subcycled timestep", false);

  // restart control
  restart_ = coordinator_list_->isParameter("restart from checkpoint file");
  if (restart_) restart_filename_ = coordinator_list_->get<std::string>("restart from checkpoint file");
}



// -----------------------------------------------------------------------------
// acquire the chosen timestep size
// -----------------------------------------------------------------------------
double Coordinator::get_dt(bool after_fail) {
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
  if (dt > max_dt_) {
    dt = max_dt_;
  }

  // ask the step manager if this step is ok
  dt = tsm_->TimeStep(S_next_->time(), dt, after_fail);
  if (subcycled_ts_) dt = std::min(dt, dt_pk);
  return dt;
}


bool Coordinator::advance(double t_old, double t_new, double& dt_next) {
  double dt = t_new - t_old;

  S_next_->advance_time(dt);
  bool fail = pk_->AdvanceStep(t_old, t_new, false);
  fail |= !pk_->ValidStep();

  // advance the iteration count and timestep size
  S_next_->advance_cycle();

  if (!fail) {
    // commit the state
    pk_->CommitStep(t_old, t_new, S_next_);

    // make observations, vis, and checkpoints
    for (const auto& obs : observations_) obs->MakeObservations(S_next_.ptr());
    visualize();
    dt_next = get_dt(fail);
    checkpoint(dt_next); // checkpoint with the new dt

    // we're done with this time step, copy the state
    *S_ = *S_next_;
    *S_inter_ = *S_next_;

  } else {
    // Failed the timestep.
    // Potentially write out failed timestep for debugging
    for (const auto& vis : failed_visualization_) WriteVis(*vis, *S_next_);

    // The timestep sizes have been updated, so copy back old soln and try again.
    *S_next_ = *S_;
    *S_inter_ = *S_;

    // check whether meshes are deformable, and if so, recover the old coordinates
    for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
         mesh!=S_->mesh_end(); ++mesh) {
      if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
        // collect the old coordinates
        std::string node_key;
        if (mesh->first != "domain") node_key= mesh->first+std::string("-vertex_coordinate");
        else node_key = std::string("vertex_coordinate");

        Teuchos::RCP<const Amanzi::CompositeVector> vc_vec = S_->GetFieldData(node_key);
        vc_vec->ScatterMasterToGhosted();
        const Epetra_MultiVector& vc = *vc_vec->ViewComponent("node", true);

        std::vector<int> node_ids(vc.MyLength());
        Amanzi::AmanziGeometry::Point_List old_positions(vc.MyLength());
        for (int n=0; n!=vc.MyLength(); ++n) {
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
    dt_next = get_dt(fail);
  }
  return fail;
}

void Coordinator::visualize(bool force) {
  // write visualization if requested
  bool dump = force;
  if (!dump) {
    for (const auto& vis : visualization_) {
      if (vis->DumpRequested(S_next_->cycle(), S_next_->time())) {
        dump = true;
      }
    }
  }

  if (dump) {
    pk_->CalculateDiagnostics(S_next_);
  }

  for (const auto& vis : visualization_) {
    if (force || vis->DumpRequested(S_next_->cycle(), S_next_->time())) {
      WriteVis(*vis, *S_next_);
    }
  }
}

void Coordinator::checkpoint(double dt, bool force) {
  if (force || checkpoint_->DumpRequested(S_next_->cycle(), S_next_->time())) {
    checkpoint_->Write(*S_next_, dt);
  }
}


// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void Coordinator::cycle_driver() {
  // wallclock duration -- in seconds
  const double duration(duration_ * 3600);

  // start at time t = t0 and initialize the state.
  double dt_restart = -1;
  {
    Teuchos::TimeMonitor monitor(*setup_timer_);
    setup();
    dt_restart = initialize();
  }

  // get the intial timestep
  double dt = dt_restart > 0 ? dt_restart : get_dt(false);

  // visualization at IC
  visualize();
  checkpoint(dt);

  // iterate process kernels
  {
    Teuchos::TimeMonitor cycle_monitor(*cycle_timer_);
#if !DEBUG_MODE
  try {
#endif
    bool fail = false;
    while ((S_->time() < t1_) &&
           ((cycle1_ == -1) || (S_->cycle() <= cycle1_)) &&
           (duration_ < 0 || timer_->totalElapsedTime(true) < duration) &&
           dt > 0.) {
      if (vo_->os_OK(Teuchos::VERB_LOW)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "======================================================================"
                  << std::endl << std::endl;
        *vo_->os() << "Cycle = " << S_->cycle();
        *vo_->os() << ",  Time [days] = "<< std::setprecision(16) << S_->time() / (60*60*24);
        *vo_->os() << ",  dt [days] = " << std::setprecision(16) << dt / (60*60*24)  << std::endl;
        *vo_->os() << "----------------------------------------------------------------------"
                  << std::endl;
      }

      *S_->GetScalarData("dt", "coordinator") = dt;
      *S_inter_->GetScalarData("dt", "coordinator") = dt;
      *S_next_->GetScalarData("dt", "coordinator") = dt;

      S_->set_initial_time(S_->time());
      S_->set_final_time(S_->time() + dt);
      S_->set_intermediate_time(S_->time());

      fail = advance(S_->time(), S_->time() + dt, dt);
    } // while not finished

#if !DEBUG_MODE
  } catch (Amanzi::Exceptions::Amanzi_exception &e) {
    // write one more vis for help debugging
    S_next_->advance_cycle();
    visualize(true); // force vis

    // flush observations to make sure they are saved
    for (const auto& obs : observations_) obs->Flush();

    // catch errors to dump two checkpoints -- one as a "last good" checkpoint
    // and one as a "debugging data" checkpoint.
    checkpoint_->set_filebasename("last_good_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), *S_, dt);
    checkpoint_->set_filebasename("error_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), *S_next_, dt);
    throw e;
  }
#endif
  }

  // finalizing simulation
  WriteStateStatistics(*S_, *vo_);
  report_memory();
  Teuchos::TimeMonitor::summarize(*vo_->os());

  finalize();
} // cycle driver


} // close namespace Amanzi
