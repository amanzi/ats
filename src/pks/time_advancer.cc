/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include <iomanip>

#include "errors.hh"
#include "Teuchos_TimeMonitor.hpp"

#include "CompositeVector.hh"
#include "TimeStepManager.hh"
#include "Visualization.hh"
#include "VisualizationDomainSet.hh"
#include "Checkpoint.hh"
#include "UnstructuredObservations.hh"
#include "State.hh"
#include "PK.hh"
#include "PK_Helpers.hh"
#include "IO.hh"

#include "time_advancer.hh"
#include "EvaluatorTimeAccumulated.hh"

namespace ATS {

TimeAdvancer::TimeAdvancer(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                           const Teuchos::RCP<Amanzi::State>& S,
                           const Teuchos::RCP<Amanzi::PK>& pk,
                           const Teuchos::RCP<Amanzi::Utils::TimeStepManager>& tsm,
                           const Amanzi::Tag& tag_current,
                           const Amanzi::Tag& tag_next,
                           const Teuchos::RCP<Amanzi::VerboseObject>& vo,
                           const Teuchos::RCP<Teuchos::Time>& wallclock_timer)
  : plist_(plist),
    vo_(vo),
    S_(S),
    pk_(pk),
    tsm_(tsm),
    tag_current_(tag_current),
    tag_next_(tag_next),
    wallclock_timer_(wallclock_timer),
    max_dt_(plist->get<double>("max timestep size [s]", 1.0e99)),
    min_dt_(plist->get<double>("min timestep size [s]", 1.0e-12)),
    cycle1_(plist->get<int>("end cycle", -1)),
    duration_(plist->get<double>("wallclock duration [hrs]", -1.0)),
    subcycled_ts_(plist->get<bool>("subcycled timestep", false))
{
  // construct checkpoint — always created so finalize() can write a final checkpoint;
  // only register with TSM for periodic checkpoints if the sublist is present.
  checkpoint_obj_ = Teuchos::rcp(
    new Amanzi::Checkpoint(plist_->sublist("checkpoint"), *S_));
  if (plist_->isSublist("checkpoint"))
    checkpoint_obj_->RegisterWithTimeStepManager(*tsm_);

  // construct observations
  if (plist_->isSublist("observations")) {
    auto& obs_list = plist_->sublist("observations");
    for (auto& entry : obs_list) {
      if (obs_list.isSublist(entry.first)) {
        observations_.emplace_back(
          Teuchos::rcp(new Amanzi::UnstructuredObservations(obs_list.sublist(entry.first))));
        observations_.back()->RegisterWithTimeStepManager(*tsm_);
      } else {
        Errors::Message msg("\"observations\" list must only include sublists.");
        Exceptions::amanzi_throw(msg);
      }
    }
  }

  // construct visualization
  if (plist_->isSublist("visualization")) {
    auto vis_list = Teuchos::sublist(plist_, "visualization");
    for (auto& entry : *vis_list) {
      std::string domain_name = entry.first;
      if (S_->HasMesh(domain_name)) {
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
        auto vis = Teuchos::rcp(new Amanzi::Visualization(*sublist_p));
        vis->set_name(domain_name);
        vis->set_mesh(mesh_p);
        vis->CreateFiles(false);
        visualization_.push_back(vis);
      } else if (Amanzi::Keys::isDomainSet(domain_name)) {
        const auto& dset = S_->GetDomainSet(Amanzi::Keys::getDomainSetName(domain_name));
        auto sublist_p = Teuchos::sublist(vis_list, domain_name);
        if (sublist_p->get("visualize individually", false)) {
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
          auto domain_name_base = Amanzi::Keys::getDomainSetName(domain_name);
          if (!sublist_p->isParameter("file name base"))
            sublist_p->set("file name base", std::string("ats_vis_") + domain_name_base);
          auto vis = Teuchos::rcp(new Amanzi::VisualizationDomainSet(*sublist_p));
          vis->set_name(domain_name_base);
          vis->set_domain_set(dset);
          vis->set_mesh(dset->getReferencingParent());
          vis->CreateFiles(false);
          visualization_.push_back(vis);
        }
      }
      for (const auto& vis : visualization_) vis->RegisterWithTimeStepManager(*tsm_);
    }
  }

  // construct failed visualization (no TSM registration needed)
  if (plist_->isSublist("visualization failed")) {
    auto vis_list = Teuchos::sublist(plist_, "visualization failed");
    for (auto& entry : *vis_list) {
      std::string domain_name = entry.first;
      if (S_->HasMesh(domain_name)) {
        auto mesh_p = S_->GetMesh(domain_name);
        auto sublist_p = Teuchos::sublist(vis_list, domain_name);
        if (!sublist_p->isParameter("file name base"))
          sublist_p->set<std::string>("file name base", std::string("ats_vis_failed_") + domain_name);
        auto vis = Teuchos::rcp(new Amanzi::Visualization(*sublist_p));
        vis->set_name(domain_name);
        vis->set_mesh(mesh_p);
        vis->CreateFiles(false);
        failed_visualization_.push_back(vis);
      }
    }
  }
}


void
TimeAdvancer::setup()
{
  S_->require_cycle(tag_next_);
  S_->Require<double>("dt", tag_next_, "dt");

  for (auto& obs : observations_) obs->Setup(S_.ptr());

  // parse "time accumulators" list: each entry is "key@tag"
  if (plist_->isParameter("time accumulators")) {
    auto keytag_strs = plist_->get<Teuchos::Array<std::string>>("time accumulators");
    for (const auto& s : keytag_strs)
      accumulator_keytags_.push_back(Amanzi::Keys::splitKeyTag(s));
  }
}


void
TimeAdvancer::initialize()
{
  // initialize dt if not already done (e.g. by a restart)
  if (!S_->GetRecord("dt", tag_next_).initialized())
    S_->GetRecordW("dt", tag_next_, "dt").set_initialized();

  // resolve time accumulator evaluators now that State is fully set up
  for (const auto& kt : accumulator_keytags_) {
    auto* eval = dynamic_cast<Amanzi::Relations::EvaluatorTimeAccumulated*>(
      &S_->GetEvaluator(kt.first, kt.second));
    if (!eval) {
      Errors::Message msg;
      msg << "TimeAdvancer: \"time accumulators\" entry \"" << kt.first << "@" << kt.second.get()
          << "\" is not an EvaluatorTimeAccumulated";
      Exceptions::amanzi_throw(msg);
    }
    accumulators_.push_back(eval);
  }

  // dump initial conditions
  visualize_();
  observe_();
  checkpoint_();
}


void
TimeAdvancer::finalize(bool checkpoint)
{
  if (checkpoint) checkpoint_obj_->Write(*S_, Amanzi::Checkpoint::WriteType::FINAL);
  for (const auto& obs : observations_) obs->Flush();
}


bool
TimeAdvancer::advance(double t_start, double t_end)
{
  // register end time as a required TSM event
  tsm_->RegisterTimeEvent(t_end);
  S_->set_time(tag_current_, t_start);
  S_->set_time(tag_next_, t_start);

  // reset all time accumulators at the start of each outer timestep
  for (auto* acc : accumulators_) acc->Reset(*S_);

  double dt = pk_->get_dt();
  if (dt > max_dt_) dt = max_dt_;
  bool fail = false;

  while (true) {
    double t_now = S_->get_time(tag_current_);

    // done?
    if (std::abs(t_end - t_now) < 1.e-10 * std::abs(t_end + 1.)) break;
    if (extraDoneCheck_(t_now, S_->get_cycle())) break;
    if (dt <= 0.) break;

    // constrain dt via TSM, then clamp
    dt = tsm_->TimeStep(t_now, dt, fail);
    if (dt < min_dt_) {
      Errors::Message msg("TimeAdvancer: timestep too small");
      Exceptions::amanzi_throw(msg);
    }
    double dt_pk = dt;
    if (subcycled_ts_) {
      double raw = pk_->get_dt();
      dt = std::min(dt, raw);
    }

    if (vo_->os_OK(Teuchos::VERB_LOW)) {
      *vo_->os() << "================================================================================"
                 << std::endl << std::endl
                 << "Cycle = " << S_->get_cycle(tag_next_)
                 << ",  Time [days] = " << std::setprecision(16)
                 << t_now / (60 * 60 * 24)
                 << ",  dt [days] = " << std::setprecision(16)
                 << dt / (60 * 60 * 24) << std::endl
                 << "--------------------------------------------------------------------------------"
                 << std::endl;
    }

    S_->Assign("dt", tag_next_, "dt", dt);
    S_->set_time(tag_next_, t_now + dt);

    fail = pk_->AdvanceStep(t_now, t_now + dt, false);

    WriteStateStatistics(*S_, *vo_, Teuchos::VERB_EXTREME);

    if (fail) {
      for (const auto& vis : failed_visualization_) WriteVis(*vis, *S_);
      pk_->FailStep(t_now, t_now + dt, tag_next_);
      S_->set_time(tag_next_, t_now);
      FailStep_(t_now, t_now + dt);
      dt = pk_->get_dt();
    } else {
      pk_->CommitStep(t_now, t_now + dt, tag_next_);
      S_->set_time(tag_current_, t_now + dt);
      S_->advance_cycle(tag_next_);
      CommitStep_(t_now, t_now + dt);
      dt = pk_->get_dt();
    }
    if (dt > max_dt_) dt = max_dt_;

    WriteStateStatistics(*S_, *vo_, Teuchos::VERB_EXTREME);
  }
  return false;
}


void
TimeAdvancer::CommitStep_(double t_old, double t_new)
{
  for (auto* acc : accumulators_) acc->Update(*S_, "time_advancer");
  visualize_();
  observe_();
  checkpoint_();
}


void
TimeAdvancer::FailStep_(double t_old, double t_new)
{
  // Deformable mesh coordinate recovery on failure is handled by the
  // VolumetricDeformation PK in its FailStep() method.
}


bool
TimeAdvancer::extraDoneCheck_(double t, int cycle) const
{
  if (cycle1_ >= 0 && cycle > cycle1_) return true;
  if (duration_ > 0. && wallclock_timer_ != Teuchos::null) {
    if (wallclock_timer_->totalElapsedTime(true) >= duration_ * 3600.) return true;
  }
  return false;
}


bool
TimeAdvancer::visualize_(bool force)
{
  bool dump = force;
  int cycle = S_->get_cycle();
  double time = S_->get_time();
  if (!dump) {
    for (const auto& vis : visualization_) dump |= vis->DumpRequested(cycle, time);
  }
  if (dump) pk_->CalculateDiagnostics(tag_next_);
  for (const auto& vis : visualization_) {
    if (force || vis->DumpRequested(cycle, time)) WriteVis(*vis, *S_);
  }
  return dump;
}


void
TimeAdvancer::observe_()
{
  for (const auto& obs : observations_) obs->MakeObservations(S_.ptr());
}


bool
TimeAdvancer::checkpoint_(bool force)
{
  int cycle = S_->get_cycle();
  double time = S_->get_time();
  bool dump = force || checkpoint_obj_->DumpRequested(cycle, time);
  if (dump) checkpoint_obj_->Write(*S_);
  return dump;
}

} // namespace ATS
