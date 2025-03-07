/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  MPC for subcycling one PK relative to another.
*/

#include "Event.hh"
#include "TimeStepManager.hh"
#include "mpc_subcycled.hh"
#include "PK_Helpers.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCSubcycled::MPCSubcycled(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_list, S, soln), MPC<PK>(pk_tree, global_list, S, soln)
{
  init_();

  // NOTE: this must be done prior to set_tags(), therefore before even parseParameterList()
  subcycling_ = plist_->get<Teuchos::Array<int>>("subcycle");
  if (subcycling_.size() != sub_pks_.size()) {
    Errors::Message msg(
      "MPCSubcycling pass \"subcycle\" list of length inconsistent with the number of PKs.\"");
    Exceptions::amanzi_throw(msg);
  }
  dts_.resize(sub_pks_.size(), -1);
  tsms_.resize(sub_pks_.size(), Teuchos::null);
}

void
MPCSubcycled::parseParameterList()
{
  for (int i = 0; i != tsms_.size(); ++i) {
    if (subcycling_[i]) {
      Teuchos::ParameterList tsm_plist(std::string("TSM: ") + sub_pks_[i]->name());
      tsms_[i] = Teuchos::rcp(new Utils::TimeStepManager(tsm_plist));

      const auto& tag = tags_[i];
      S_->require_time(tag.first);
      S_->require_time(tag.second);
    }
  }

  MPC<PK>::parseParameterList();

  // min dt allowed in subcycling
  target_dt_ = plist_->get<double>("subcycling target timestep [s]", -1);
}


void
MPCSubcycled::set_tags(const Tag& current, const Tag& next)
{
  PK::set_tags(current, next);

  tags_.clear();
  int i = 0;
  for (auto& pk : sub_pks_) {
    if (subcycling_[i]) {
      Tag lcurrent(Keys::cleanName(pk->name() + " current"));
      Tag lnext(Keys::cleanName(pk->name() + " next"));
      tags_.emplace_back(std::make_pair(lcurrent, lnext));
    } else {
      tags_.emplace_back(std::make_pair(current, next));
    }
    pk->set_tags(tags_.back().first, tags_.back().second);
    ++i;
  }
}


void
MPCSubcycled::Setup()
{
  int i = 0;
  for (const auto& tag : tags_) {
    S_->require_cycle(tag.second);
    if (subcycling_[i]) S_->Require<double>("dt", tag.first, name());
    ++i;
  }
  MPC<PK>::Setup();
}

void
MPCSubcycled::Initialize()

{
  int i = 0;
  for (const auto& tag : tags_) {
    if (subcycling_[i]) S_->GetRecordW("dt", tag.first, name()).set_initialized();
    ++i;
  }
  MPC<PK>::Initialize();
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
MPCSubcycled::get_dt()
{
  double dt = std::numeric_limits<double>::max();
  if (target_dt_ > 0) dt = target_dt_;

  int i = 0;
  for (auto& pk : sub_pks_) {
    dts_[i] = pk->get_dt();
    if (!subcycling_[i]) dt = std::min(dt, dts_[i]);
    ++i;
  }

  dt_ = dt;
  return dt;
}


// -----------------------------------------------------------------------------
// Set standard dt
// -----------------------------------------------------------------------------
void
MPCSubcycled::set_dt(double dt)
{
  dt_ = dt;
  int i = 0;
  for (auto& pk : sub_pks_) {
    if (!subcycling_[i] || dt < dts_[i]) {
      dts_[i] = dt;
      pk->set_dt(dt);
    }
    ++i;
  }
}


// -----------------------------------------------------------------------------
// Advance the ith sub-PK individually, returning a failure
// -----------------------------------------------------------------------------
bool
MPCSubcycled::AdvanceStep_i_(std::size_t i, double t_old, double t_new, bool reinit)
{
  bool fail = false;
  if (subcycling_[i]) {
    // advance the subcycled, subcycling if needed
    bool done = false;
    double t_inner = t_old;
    double dt_inner = dts_[i];
    bool fail_inner = false;

    Tag tag_subcycle_current = tags_[i].first;
    Tag tag_subcycle_next = tags_[i].second;

    tsms_[i]->RegisterTimeEvent(t_new);
    S_->set_time(tag_subcycle_current, t_old);

    while (!done) {
      dt_inner = tsms_[i]->TimeStep(t_inner, dt_inner, fail_inner);
      S_->Assign("dt", tag_subcycle_current, name(), dt_inner);
      S_->set_time(tag_subcycle_next, t_inner + dt_inner);
      fail_inner = sub_pks_[i]->AdvanceStep(t_inner, t_inner + dt_inner, false);

      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  step failed? " << fail_inner << std::endl;

      if (fail_inner) {
        sub_pks_[i]->FailStep(t_old, t_new, tag_subcycle_next);
        dt_inner = sub_pks_[i]->get_dt();

        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  failed, new timestep is " << dt_inner << std::endl;

      } else {
        sub_pks_[i]->CommitStep(t_inner, t_inner + dt_inner, tag_subcycle_next);
        t_inner += dt_inner;
        if (std::abs(t_new - t_inner) < 1.e-10) done = true;

        S_->set_time(tag_subcycle_current, S_->get_time(tag_subcycle_next));
        S_->advance_cycle(tag_subcycle_next);

        dt_inner = sub_pks_[i]->get_dt();
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  success, new timestep is " << dt_inner << std::endl;
      }
    }

  } else {
    // advance the standard PK using the full step size
    fail = sub_pks_[i]->AdvanceStep(t_old, t_new, reinit);
  }
  return fail;
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
MPCSubcycled::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;

  int i = 0;
  for (auto& pk : sub_pks_) {
    fail = AdvanceStep_i_(i, t_old, t_new, reinit);
    if (fail) return fail;
    ++i;
  }

  return fail;
}


void
MPCSubcycled::CommitStep(double t_old, double t_new, const Tag& tag)
{
  MPC<PK>::CommitStep(t_old, t_new, tag);

  if (S_->get_cycle() < 0 && tag == Tags::NEXT) {
    // initial commit, also do the substep commits
    int i = 0;
    for (auto& pk : sub_pks_) {
      if (subcycling_[i]) {
        pk->CommitStep(t_old, t_new, tags_[i].second);
      }
      ++i;
    }
  }
}


} // namespace Amanzi
