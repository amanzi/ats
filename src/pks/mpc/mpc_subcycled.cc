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

#include "TimeStepManager.hh"
#include "time_advancer.hh"
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
  time_advancers_.resize(sub_pks_.size(), Teuchos::null);
}


void
MPCSubcycled::parseParameterList()
{
  for (int i = 0; i != time_advancers_.size(); ++i) {
    if (subcycling_[i]) {
      // Locate the TimeAdvancer sublist: "PKNAME time advancer" if present,
      // else copy from shared "time advancer" if present, else create empty.
      // Using Teuchos::sublist ensures the sublist is created in plist_ so the
      // final parameters used are recorded there.
      std::string pk_ta_key = sub_pks_[i]->name() + " time advancer";
      if (!plist_->isSublist(pk_ta_key) && plist_->isSublist("time advancer"))
        plist_->sublist(pk_ta_key) = plist_->sublist("time advancer");
      auto ta_plist = Teuchos::sublist(plist_, pk_ta_key);

      // create a TSM for this subcycled PK
      auto tsm = Teuchos::rcp(new Utils::TimeStepManager());

      const auto& tag = tags_[i];
      S_->require_time(tag.first);
      S_->require_time(tag.second);

      time_advancers_[i] = Teuchos::rcp(new ATS::TimeAdvancer(
        ta_plist, S_, sub_pks_[i], tsm,
        tag.first, tag.second, vo_));
    }
  }

  MPC<PK>::parseParameterList();

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
    ++i;
  }
  MPC<PK>::Setup();
  for (auto& ta : time_advancers_) {
    if (ta.get()) ta->setup();
  }
}


void
MPCSubcycled::Initialize()
{
  int i = 0;
  for (const auto& tag : tags_) {
    if (subcycling_[i]) {
      S_->set_time(tag.first, S_->get_time(Amanzi::Tags::CURRENT));
      S_->set_time(tag.second, S_->get_time(Amanzi::Tags::NEXT));
    }
    ++i;
  }
  MPC<PK>::Initialize();
  for (auto& ta : time_advancers_) {
    if (ta.get()) ta->initialize();
  }
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
  if (subcycling_[i]) {
    return time_advancers_[i]->advance(t_old, t_new);
  } else {
    return sub_pks_[i]->AdvanceStep(t_old, t_new, reinit);
  }
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
