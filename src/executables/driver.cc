/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include <iomanip>
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
#include "MeshInfo.hh"
#include "IO.hh"
#include "GeometricModel.hh"
#include "Checkpoint.hh"
#include "State.hh"
#include "PK.hh"
#include "TreeVector.hh"
#include "PK_Factory.hh"
#include "PK_Helpers.hh"

#include "Key.hh"
#include "ats_mesh_factory.hh"
#include "time_advancer.hh"
#include "driver.hh"

namespace ATS {

Driver::Driver(const Teuchos::RCP<Teuchos::ParameterList>& plist,
               const Teuchos::RCP<Teuchos::Time>& wallclock_timer,
               const Teuchos::RCP<const Teuchos::Comm<int>>& teuchos_comm,
               const Amanzi::Comm_ptr_type& comm)
  : plist_(plist),
    teuchos_comm_(teuchos_comm),
    comm_(comm),
    restart_(false),
    wallclock_timer_(wallclock_timer),
    wallclock_monitor_(*wallclock_timer)
{
  coordinator_list_ = Teuchos::sublist(plist_, "cycle driver");
  vo_ = Teuchos::rcp(new Amanzi::VerboseObject(comm_, "ATS", *coordinator_list_));
  Teuchos::OSTab tab = vo_->getOSTab();

  timers_["0: create mesh"] = Teuchos::TimeMonitor::getNewCounter("0: create mesh");
  timers_["1: create run"] = Teuchos::TimeMonitor::getNewCounter("1: create run");
  timers_["2a: parseParameterList"] = Teuchos::TimeMonitor::getNewCounter("2a: parseParameterList");
  timers_["2b: setup"] = Teuchos::TimeMonitor::getNewCounter("2b: setup");
  timers_["3: initialize"] = Teuchos::TimeMonitor::getNewCounter("3: initialize");
  timers_["4: solve"] = Teuchos::TimeMonitor::getNewCounter("4: solve");
  timers_["4a: advance step"] = Teuchos::TimeMonitor::getNewCounter("4a: advance step");
  timers_["4b: visualize"] = Teuchos::TimeMonitor::getNewCounter("4b: visualize");
  timers_["4c: observe"] = Teuchos::TimeMonitor::getNewCounter("4c: observe");
  timers_["4d: checkpoint"] = Teuchos::TimeMonitor::getNewCounter("4d: checkpoint");
  timers_["5: finalize"] = Teuchos::TimeMonitor::getNewCounter("5: finalize");

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
    S_ = Teuchos::rcp(new Amanzi::State(plist_->sublist("state")));

    Teuchos::ParameterList reg_list = plist_->sublist("regions");
    Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
      Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_list, *comm_));

    ATS::Mesh::createMeshes(plist_, comm_, gm, *S_);
  }
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("0: create mesh");
  }

  // create PKs and supporting objects
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "================================================================================"
               << std::endl
               << "Beginning run creation stage..." << std::endl
               << std::flush;
  }
  {
    Teuchos::TimeMonitor timer(*timers_.at("1: create run"));
    initializeFromPlist_();

    Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(plist_, "PKs");
    Teuchos::ParameterList pk_tree_list = coordinator_list_->sublist("PK tree");
    if (pk_tree_list.numParams() != 1) {
      Errors::Message message("Driver: PK tree list should contain exactly one root node list");
      Exceptions::amanzi_throw(message);
    }
    Teuchos::ParameterList::ConstIterator pk_item = pk_tree_list.begin();
    const std::string& pk_name = pk_tree_list.name(pk_item);

    soln_ = Teuchos::rcp(new Amanzi::TreeVector(comm_));

    Amanzi::PKFactory pk_factory;
    pk_ = pk_factory.CreatePK(pk_name, pk_tree_list, plist_, S_, soln_);

    for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
      auto& mesh_sublist = plist_->sublist("mesh").sublist(mesh->first);
      if (mesh_sublist.isSublist("region analysis")) {
        Amanzi::InputAnalysis analysis(mesh->second.first, mesh->first);
        analysis.Init(mesh_sublist.sublist("region analysis"));
        analysis.RegionAnalysis();
        analysis.OutputBCs();
      }
      if (mesh_sublist.isSublist("mesh info")) {
        Teuchos::RCP<Amanzi::MeshInfo> mesh_info =
          Teuchos::rcp(new Amanzi::MeshInfo(mesh_sublist.sublist("mesh info"), *S_));
        mesh_info->WriteMeshCentroids(mesh->first, *(mesh->second.first));
      }
    }

    // create the TimeStepManager
    tsm_ = Teuchos::rcp(
      new Amanzi::Utils::TimeStepManager(coordinator_list_->sublist("timestep manager")));
    tsm_->RegisterTimeEvent(t1_);
    if (coordinator_list_->isSublist("required times")) {
      Teuchos::ParameterList& sublist = coordinator_list_->sublist("required times");
      Amanzi::Utils::IOEvent pause_times(sublist);
      pause_times.RegisterWithTimeStepManager(*tsm_);
    }
  }
  if (vo_->os_OK(Teuchos::VERB_LOW)) {
    *vo_->os() << "  ... completed: ";
    reportOneTimer_("1: create run");
  }
}


void
Driver::parseParameterList()
{
  Teuchos::TimeMonitor monitor(*timers_.at("2a: parseParameterList"));

  S_->require_time(Amanzi::Tags::CURRENT);
  S_->require_time(Amanzi::Tags::NEXT);

  pk_->set_tags(Amanzi::Tags::CURRENT, Amanzi::Tags::NEXT);
  pk_->parseParameterList();
}


void
Driver::setup()
{
  Teuchos::TimeMonitor monitor(*timers_.at("2b: setup"));

  S_->Require<double>("atmospheric_pressure", Amanzi::Tags::DEFAULT, "coordinator");
  S_->Require<Amanzi::AmanziGeometry::Point>("gravity", Amanzi::Tags::DEFAULT, "coordinator");

  // require vertex coordinates on Tags::NEXT for all deformable meshes so that
  // Driver can initialize, suppress vis output, and handle checkpoint restart
  // independently of any PK that may also require them on a different tag.
  for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
    if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
      Amanzi::Key node_key = Amanzi::Keys::getKey(mesh->first, "vertex_coordinates");
      S_->Require<Amanzi::CompositeVector, Amanzi::CompositeVectorSpace>(
        node_key, Amanzi::Tags::NEXT, node_key)
        .SetMesh(mesh->second.first)
        ->SetComponent("node", Amanzi::AmanziMesh::Entity_kind::NODE, mesh->second.first->getSpaceDimension());
    }
  }

  // order matters: PK leaves first, then observations can use those fields,
  // then State::Setup() allocates memory for all secondaries
  pk_->Setup();
  if (time_advancer_.get()) time_advancer_->setup();
  S_->Setup();
}


void
Driver::initialize()
{
  Teuchos::TimeMonitor monitor(*timers_.at("3: initialize"));
  Teuchos::OSTab tab = vo_->getOSTab();

  S_->set_time(Amanzi::Tags::CURRENT, t0_);
  S_->set_time(Amanzi::Tags::NEXT, t0_);
  S_->set_cycle(cycle0_);

  if (restart_) {
    double t_restart = Amanzi::ReadCheckpointInitialTime(comm_, restart_filename_);
    for (const auto& record : S_->GetRecordSet("time")) {
      S_->set_time(record.first, t_restart);
    }
    t0_ = t_restart;
  }

  S_->InitializeFields();
  pk_->Initialize();
  pk_->CommitStep(t0_, t0_, Amanzi::Tags::NEXT);

  // initialize vertex coordinates and suppress vis (checkpoint handles them)
  for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
    if (S_->IsDeformableMesh(mesh->first) && !S_->IsAliasedMesh(mesh->first)) {
      Amanzi::Key node_key = Amanzi::Keys::getKey(mesh->first, "vertex_coordinates");
      copyMeshCoordinatesToVector(
        *mesh->second.first,
        S_->GetW<Amanzi::CompositeVector>(node_key, Amanzi::Tags::NEXT, node_key));
      S_->GetRecordW(node_key, Amanzi::Tags::NEXT, node_key).set_initialized();
      S_->GetRecordW(node_key, Amanzi::Tags::NEXT, node_key).set_io_vis(false);
    }
  }

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

  S_->InitializeEvaluators();
  S_->InitializeFieldCopies();
  S_->CheckAllFieldsInitialized();
  S_->InitializeIOFlags();

  pk_->CommitStep(S_->get_time(), S_->get_time(), Amanzi::Tags::NEXT);

  if (S_->get_cycle() == -1) S_->advance_cycle();

  if (time_advancer_.get()) time_advancer_->initialize();
}


void
Driver::finalize(bool checkpoint)
{
  Teuchos::TimeMonitor monitor(*timers_.at("5: finalize"));

  pk_->CalculateDiagnostics(Amanzi::Tags::NEXT);
  if (time_advancer_.get()) time_advancer_->finalize(checkpoint);
  WriteStateStatistics(*S_, *vo_);
  report_memory();
}


static double
rss_usage()
{
#if (defined(__unix__) || defined(__unix) || defined(unix) || defined(__APPLE__) || \
     defined(__MACH__))
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#if (defined(__APPLE__) || defined(__MACH__))
  return static_cast<double>(usage.ru_maxrss) / 1024.0 / 1024.0;
#else
  return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
#else
  return 0.0;
#endif
}


void
Driver::report_memory()
{
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    double global_ncells(0.0);
    double local_ncells(0.0);
    for (Amanzi::State::mesh_iterator mesh = S_->mesh_begin(); mesh != S_->mesh_end(); ++mesh) {
      Epetra_Map cell_map =
        (mesh->second.first)->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false);
      global_ncells += cell_map.NumGlobalElements();
      local_ncells += cell_map.NumMyElements();
    }

    double mem = rss_usage();
    double percell = (local_ncells > 0) ? mem / local_ncells : mem;

    double max_percell(0.0), min_percell(0.0);
    comm_->MinAll(&percell, &min_percell, 1);
    comm_->MaxAll(&percell, &max_percell, 1);

    double total_mem(0.0), max_mem(0.0), min_mem(0.0);
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
}


void
Driver::initializeFromPlist_()
{
  Amanzi::Utils::Units units;
  t0_ = coordinator_list_->get<double>("start time");
  std::string t0_units = coordinator_list_->get<std::string>("start time units", "s");
  if (!units.IsValidTime(t0_units)) {
    Errors::Message msg;
    msg << "Driver start time: unknown time units type: \"" << t0_units
        << "\"  Valid are: " << units.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
  bool success;
  t0_ = units.ConvertTime(t0_, t0_units, "s", success);

  t1_ = coordinator_list_->get<double>("end time");
  std::string t1_units = coordinator_list_->get<std::string>("end time units", "s");
  if (!units.IsValidTime(t1_units)) {
    Errors::Message msg;
    msg << "Driver end time: unknown time units type: \"" << t1_units
        << "\"  Valid are: " << units.ValidTimeStrings();
    Exceptions::amanzi_throw(msg);
  }
  t1_ = units.ConvertTime(t1_, t1_units, "s", success);

  cycle0_ = coordinator_list_->get<int>("start cycle", -1);

  restart_ = coordinator_list_->isParameter("restart from checkpoint file");
  if (restart_)
    restart_filename_ = coordinator_list_->get<std::string>("restart from checkpoint file");
}


void
Driver::reportOneTimer_(const std::string& timer_name)
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
