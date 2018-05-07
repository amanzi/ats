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

#include "InputAnalysis.hh"
#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"
//#include "pk_factory_ats.hh"

#include "coordinator.hh"

#define DEBUG_MODE 1

namespace ATS {

Coordinator::Coordinator(Teuchos::ParameterList& parameter_list,
                         Teuchos::RCP<Amanzi::State>& S,
                         Epetra_MpiComm* comm ) :
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
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(pks_list, pk_name);
  pk_list->set("PK name", pk_name);

  pk_ = pk_factory.CreatePK(pk_name, pk_tree_list, parameter_list_, S_, soln_);

  int rank = comm_->MyPID();
  int size = comm_->NumProc();
  std::stringstream check;
  
  if(parameter_list_->sublist("mesh").isSublist("column"))  
    check << "checkpoint " << rank;
  else
    check << "checkpoint";
  
  // create the checkpointing

  Teuchos::ParameterList& chkp_plist = parameter_list_->sublist(check.str());
  if (parameter_list_->sublist("mesh").isSublist("column") && size >1){
    MPI_Comm mpi_comm_self(MPI_COMM_SELF);
    Epetra_MpiComm *comm_self = new Epetra_MpiComm(mpi_comm_self);
    checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(chkp_plist, comm_self));
  }
  else
    checkpoint_ = Teuchos::rcp(new Amanzi::Checkpoint(chkp_plist, comm_));
  

  // create the observations
  Teuchos::ParameterList& observation_plist = parameter_list_->sublist("observations");
  observations_ = Teuchos::rcp(new Amanzi::UnstructuredObservations(observation_plist,
          Teuchos::null, comm_));

  // check whether meshes are deformable, and if so require a nodal position
  for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
       mesh!=S_->mesh_end(); ++mesh) {

    if (S_->IsDeformableMesh(mesh->first) ){
      if (mesh->first.find("column") != std::string::npos) {
        std::string node_key = mesh->first+std::string("-vertex_coordinate");
        S_->RequireField(node_key)->SetMesh(mesh->second.first)->SetGhosted()
          ->AddComponent("node", Amanzi::AmanziMesh::NODE, mesh->second.first->space_dimension());
      }
      else if (!parameter_list_->sublist("mesh").isSublist("column")){
        std::string node_key;
        if (mesh->first != "domain")
          node_key= mesh->first+std::string("-vertex_coordinate");
        else
          node_key = std::string("vertex_coordinate");

        S_->RequireField(node_key)->SetMesh(mesh->second.first)->SetGhosted()
          ->AddComponent("node", Amanzi::AmanziMesh::NODE, mesh->second.first->space_dimension()); 
      }
    }

    
  // -------------- ANALYSIS --------------------------------------------
    if (parameter_list_->isSublist("analysis")){
      Amanzi::InputAnalysis analysis(mesh->second.first);
      //analysis.Init(parameter_list_->sublist(mesh->first));
      //std::cout<<mesh->first<<"\n";
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

  pk_->Setup(S_.ptr());  
  S_->Setup();
}

void Coordinator::initialize() {
  // Restart from checkpoint, part 1.

  // This is crufty -- blame the BDF1 time integrator, whose solution history
  // needs to be updated to accept the new time as its initial time.

  // Note that if this is so, we can probably ignore some of the above
  // initialize() calls and the commit_state() call, but I'm afraid to try
  // that and break all the PKs.
  // Currently not a true restart -- for a true restart this should also get:
  // -- timestep size dt
  // -- BDF history to allow projection to continue correctly.

  int size = comm_->NumProc();


  //---
  if (restart_) {
    if (parameter_list_->sublist("mesh").isSublist("column") && size >1){
      MPI_Comm mpi_comm_self(MPI_COMM_SELF);
      Epetra_MpiComm *comm_self = new Epetra_MpiComm(mpi_comm_self);
      t0_ = Amanzi::ReadCheckpointInitialTime(comm_self, restart_filename_);
      S_->set_time(t0_);
    }
    else{
      t0_ = Amanzi::ReadCheckpointInitialTime(comm_, restart_filename_);
      S_->set_time(t0_);
    }
  }

  // Restart from checkpoint, part 2.
  if (restart_) {
    if (parameter_list_->sublist("mesh").isSublist("column") && size >1){
      MPI_Comm mpi_comm_self(MPI_COMM_SELF);
      Epetra_MpiComm *comm_self = new Epetra_MpiComm(mpi_comm_self);
      ReadCheckpoint(comm_self, S_.ptr(), restart_filename_);
      t0_ = S_->time();
      cycle0_ = S_->cycle();
    }
    else{
      ReadCheckpoint(comm_, S_.ptr(), restart_filename_);
      t0_ = S_->time();
      cycle0_ = S_->cycle();
    }
    
    for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
         mesh!=S_->mesh_end(); ++mesh) {
      if (boost::starts_with(mesh->first, "column")){
        DeformCheckpointMesh(S_.ptr(), mesh->first);
      }
    }
    
  }
  
  // double check columns
  if(restart_){
    for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
         mesh!=S_->mesh_end(); ++mesh) {
      if (boost::starts_with(mesh->first, "column")){
        DeformCheckpointMesh(S_.ptr(), mesh->first);
      }
    }
  }

  
  // Initialize the state (initializes all dependent variables).
  //S_->Initialize();
  *S_->GetScalarData("dt", "coordinator") = 0.;
  S_->GetField("dt","coordinator")->set_initialized();

  S_->InitializeFields();

  // Initialize the process kernels (initializes all independent variables)
  pk_->Initialize(S_.ptr());

 // Final checks.
  S_->CheckNotEvaluatedFieldsInitialized();

  S_->InitializeEvaluators();


  S_->CheckAllFieldsInitialized();

  S_->WriteStatistics(vo_);


  // commit the initial conditions.
  pk_->CommitStep(0., 0., S_);

  // visualization
  auto vis_list = Teuchos::sublist(parameter_list_,"visualization");
  for (auto& entry : *vis_list) {
    std::string domain_name = entry.first;

    if (S_->HasMesh(domain_name)) {
      // visualize standard domain
      auto mesh_p = S_->GetMesh(domain_name);
      if (domain_name == "surface" && S_->HasMesh("surface_3d"))
        mesh_p = S_->GetMesh("surface_3d");
      
      auto sublist_p = Teuchos::sublist(vis_list, domain_name);

      // vis successful timesteps
      auto vis = Teuchos::rcp(new Amanzi::Visualization(*sublist_p));
      vis->set_name(domain_name);
      vis->set_mesh(mesh_p);
      vis->CreateFiles();
    
      visualization_.push_back(vis);

    } else if (boost::ends_with(domain_name, "_*")) {
      // visualize domain set
      std::string domain_set_name = domain_name.substr(0,domain_name.size()-2);
      for (auto m=S_->mesh_begin(); m!=S_->mesh_end(); ++m) {
        if (boost::starts_with(m->first, domain_set_name)) {
          // visualize each subdomain
          Teuchos::ParameterList sublist = vis_list->sublist(domain_name);
          sublist.set<std::string>("file name base", std::string("visdump_")+m->first);
          auto vis = Teuchos::rcp(new Amanzi::Visualization(sublist));
          vis->set_name(m->first);
          vis->set_mesh(m->second.first);    
          vis->CreateFiles();
          visualization_.push_back(vis);
        }
      }

    }
  }

  // make observations
  observations_->MakeObservations(*S_);

  S_->set_time(t0_); // in case steady state solve changed this
  S_->set_cycle(cycle0_);

  // set up the TSM
  // -- register visualization times
  for (std::vector<Teuchos::RCP<Amanzi::Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    (*vis)->RegisterWithTimeStepManager(tsm_.ptr());
  }

  // -- register checkpoint times
  checkpoint_->RegisterWithTimeStepManager(tsm_.ptr());

  // -- register observation times
  observations_->RegisterWithTimeStepManager(tsm_.ptr());

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
  S_inter_ = Teuchos::rcp(new Amanzi::State(*S_));
  *S_inter_ = *S_;

  // set the states in the PKs
  //Teuchos::RCP<const State> cS = S_; // ensure PKs get const reference state
  pk_->set_states(S_, S_inter_, S_next_); // note this does not allow subcycling


  
}

void Coordinator::finalize() {
  // Force checkpoint at the end of simulation.
  // Only do if the checkpoint was not already written, or we would be writing
  // the same file twice.
  // This really should be removed, but for now is left to help stupid developers.
  if (!checkpoint_->DumpRequested(S_next_->cycle(), S_next_->time())) {
    pk_->CalculateDiagnostics(S_next_);
    WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr(), 0.0);
  }

  // flush observations to make sure they are saved
  observations_->Flush();
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
    *vo_->os() << "  Minumum per core:   " << std::setw(7) << min_mem 
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
  t0_ = coordinator_list_->get<double>("start time");
  t1_ = coordinator_list_->get<double>("end time");
  std::string t0_units = coordinator_list_->get<std::string>("start time units", "s");
  std::string t1_units = coordinator_list_->get<std::string>("end time units", "s");

  if (t0_units == "s") {
    // internal units in s
  } else if (t0_units == "d") { // days
    t0_ = t0_ * 24.0*3600.0;
  } else if (t0_units == "yr") { // years
    t0_ = t0_ * 365.*24.0*3600.0;
  } else {
    Errors::Message message("Coordinator: error, invalid start time units");
    Exceptions::amanzi_throw(message);
  }

  if (t1_units == "s") {
    // internal units in s
  } else if (t1_units == "d") { // days
    t1_ = t1_ * 24.0*3600.0;
  } else if (t1_units == "yr") { // years
    t1_ = t1_ * 365.25*24.0*3600.0;
  } else {
    Errors::Message message("Coordinator: error, invalid end time units");
    Exceptions::amanzi_throw(message);
  }

  max_dt_ = coordinator_list_->get<double>("max time step size", 1.0e99);
  min_dt_ = coordinator_list_->get<double>("min time step size", 1.0e-12);
  cycle0_ = coordinator_list_->get<int>("start cycle",0);
  cycle1_ = coordinator_list_->get<int>("end cycle",-1);
  duration_ = coordinator_list_->get<double>("wallclock duration [hrs]", -1.0);

  // restart control
  restart_ = coordinator_list_->isParameter("restart from checkpoint file");
  if (restart_) {
    restart_filename_ = coordinator_list_->get<std::string>("restart from checkpoint file");
    // likely should ensure the file exists here? --etc
  }

  int rank = comm_->MyPID();
  
  if (coordinator_list_->isSublist("restart from checkpoint columns")){
    restart_ = true;
    Teuchos::ParameterList list = coordinator_list_->sublist("restart from checkpoint columns");
    std::stringstream name;
    name << list.get<std::string>("checkpoint file") << "checkpoint_" << rank << "_" << list.get<std::string>("cycles") << ".h5";
    restart_filename_  = name.str();
  }
}



// -----------------------------------------------------------------------------
// acquire the chosen timestep size
// -----------------------------------------------------------------------------
double Coordinator::get_dt(bool after_fail) {
  // get the physical step size
  double dt = pk_->get_dt();

  if (dt < 0.) {
    return dt;
  }

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
  return dt;
}


bool Coordinator::advance(double t_old, double t_new) {
  double dt = t_new - t_old;

  S_next_->advance_time(dt);
  bool fail = pk_->AdvanceStep(t_old, t_new, false);
  fail |= !pk_->ValidStep();

  // advance the iteration count and timestep size
  S_next_->advance_cycle();

  if (!fail) {
    // commit the state

    if (coordinator_list_->get<bool>("subcycle", false)){}
    else{
        pk_->CommitStep(t_old, t_new, S_next_);
      }

    // make observations, vis, and checkpoints
    observations_->MakeObservations(*S_next_);
    visualize();
    checkpoint(dt);

    // we're done with this time step, copy the state
    *S_ = *S_next_;
    *S_inter_ = *S_next_;

  } else {
    // Failed the timestep.  
    // Potentially write out failed timestep for debugging
    for (std::vector<Teuchos::RCP<Amanzi::Visualization> >::iterator vis=failed_visualization_.begin();
         vis!=failed_visualization_.end(); ++vis) {
      WriteVis((*vis).ptr(), S_next_.ptr());
    }

    // The timestep sizes have been updated, so copy back old soln and try again.
    *S_next_ = *S_;

    // check whether meshes are deformable, and if so, recover the old coordinates
    for (Amanzi::State::mesh_iterator mesh=S_->mesh_begin();
         mesh!=S_->mesh_end(); ++mesh) {
      bool surf = boost::starts_with(mesh->first, "surface_");

      if (S_->IsDeformableMesh(mesh->first)){
        if (mesh->first.find("column") != std::string::npos) {
        // collect the old coordinates
          
          std::string node_key = mesh->first+std::string("-vertex_coordinate");
          
          Teuchos::RCP<const Amanzi::CompositeVector> vc_vec = S_->GetFieldData(node_key);
          vc_vec->ScatterMasterToGhosted();
          const Epetra_MultiVector& vc = *vc_vec->ViewComponent("node", true);
          std::vector<int> node_ids(vc.MyLength());
          Amanzi::AmanziGeometry::Point_List old_positions(vc.MyLength());
          for (int n=0;n!=vc.MyLength();++n) {
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
        
        else if (!parameter_list_->sublist("mesh").isSublist("column")) {
          // collect the old coordinates
          
          std::string node_key;
          if (mesh->first != "domain")
            node_key= mesh->first+std::string("-vertex_coordinate");
          else
            node_key = std::string("vertex_coordinate");

          Teuchos::RCP<const Amanzi::CompositeVector> vc_vec = S_->GetFieldData(node_key);
          vc_vec->ScatterMasterToGhosted();
          const Epetra_MultiVector& vc = *vc_vec->ViewComponent("node", true);
          std::vector<int> node_ids(vc.MyLength());
          Amanzi::AmanziGeometry::Point_List old_positions(vc.MyLength());
          for (int n=0;n!=vc.MyLength();++n) {
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
  }
  return fail;
}

void Coordinator::visualize(bool force) {
  // write visualization if requested
  bool dump = force;
  if (!dump) {
    for (std::vector<Teuchos::RCP<Amanzi::Visualization> >::iterator vis=visualization_.begin();
         vis!=visualization_.end(); ++vis) {
      if ((*vis)->DumpRequested(S_next_->cycle(), S_next_->time())) {
        dump = true;
      }
    }
  }

  if (dump) {
    pk_->CalculateDiagnostics(S_next_);
  }

  for (std::vector<Teuchos::RCP<Amanzi::Visualization> >::iterator vis=visualization_.begin();
       vis!=visualization_.end(); ++vis) {
    if (force || (*vis)->DumpRequested(S_next_->cycle(), S_next_->time())) {
      WriteVis((*vis).ptr(), S_next_.ptr());
    }
  }
}

void Coordinator::checkpoint(double dt, bool force) {
  if (force || checkpoint_->DumpRequested(S_next_->cycle(), S_next_->time())) {
    WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr(), dt);
  }
}


// -----------------------------------------------------------------------------
// timestep loop
// -----------------------------------------------------------------------------
void Coordinator::cycle_driver() {
  // wallclock duration -- in seconds
  const double duration(duration_ * 3600);

  // start at time t = t0 and initialize the state.
  {
    Teuchos::TimeMonitor monitor(*setup_timer_);
    setup();   
    initialize();
   
  }

  //  exit(0);

  // get the intial timestep -- note, this would have to be fixed for a true restart
  double dt = get_dt(false);

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
      if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *vo_->os() << "======================================================================"
                  << std::endl << std::endl;
        *vo_->os() << "Cycle = " << S_->cycle();
        *vo_->os() << std::setprecision(15) << ",  Time [days] = "<< S_->time() / (60*60*24);
        *vo_->os() << ",  dt [days] = " << dt / (60*60*24)  << std::endl;
        *vo_->os() << "----------------------------------------------------------------------"
                  << std::endl;
      }

      *S_->GetScalarData("dt", "coordinator") = dt;
      *S_inter_->GetScalarData("dt", "coordinator") = dt;
      *S_next_->GetScalarData("dt", "coordinator") = dt;

      S_->set_initial_time(S_->time());
      S_->set_final_time(S_->time() + dt);
      S_->set_intermediate_time(S_->time());

      fail = advance(S_->time(), S_->time() + dt);
      //S_->WriteStatistics(vo_);  
      dt = get_dt(fail);

    } // while not finished


#if !DEBUG_MODE
  }

  catch (Amanzi::Exceptions::Amanzi_exception &e) {
    // write one more vis for help debugging
    S_next_->advance_cycle();
    visualize(true); // force vis

    // flush observations to make sure they are saved
    observations_->Flush();

    // catch errors to dump two checkpoints -- one as a "last good" checkpoint
    // and one as a "debugging data" checkpoint.
    checkpoint_->set_filebasename("last_good_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), S_.ptr(), dt);
    checkpoint_->set_filebasename("error_checkpoint");
    WriteCheckpoint(checkpoint_.ptr(), S_next_.ptr(), dt);
    throw e;
  }
#endif
  }


  // finalizing simulation                                                                                                                                                                                                               
  S_->WriteStatistics(vo_);  
  report_memory();
  Teuchos::TimeMonitor::summarize(*vo_->os());

  finalize();

} // cycle driver

} // close namespace Amanzi
