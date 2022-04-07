/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  MPC for subcycling one PK relative to another.

  NOTE: this is currently a hack-job, as it really does flow + transport
  coupling.  It will be made more general, and an MPC flow + transport will be
  done better eventually, but for now, we proceed forward.  Blocked by ATS#115.
*/

#include "mpc_subcycled.hh"
#include "pk_helpers.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCSubcycled::MPCSubcycled(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln) :
  PK(pk_tree, global_list, S, soln),
  MPC<PK>(pk_tree, global_list, S, soln)
{
  init_();

  // Master PK is the PK whose time step size sets the size, the subcycled is subcycled.
  subcycling_ = plist_->get<Teuchos::Array<int>>("subcycle");
  if (subcycling_.size() != sub_pks_.size()) {
    Errors::Message msg("MPCSubcycling pass \"subcycle\" list of length inconsistent with the number of PKs.\"");
    Exceptions::amanzi_throw(msg);
  }
  dts_.resize(sub_pks_.size(), -1);

  // min dt allowed in subcycling
  min_dt_ = plist_->get<double>("minimum subcycled relative dt", 1.e-5);
}


void
MPCSubcycled::set_tags(const Tag& current, const Tag& next)
{
  PK::set_tags(current, next);

  tags_.clear();
  int i = 0;
  for (auto& pk : sub_pks_) {
    if (subcycling_[i]) {
      tags_.emplace_back(std::make_pair(Tag{pk->name()+"_current"}, Tag{pk->name()+"_next"}));
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
    S_->require_time(tag.first);
    S_->require_time(tag.second);
    S_->require_cycle(tag.second);
    if (subcycling_[i]) S_->Require<double>("dt", tag.second, name());
    ++i;
  }
  sub_pks_[0]->Setup();

  // BEGIN HACK --etc
  // hack -- assign water flux eval and field to transport's next
  // This mimics the concept of "pointer" evaluators and fields.
  // aliasVector(*S_, "water_flux", tag_next_, tags_[1].second);
  // aliasVector(*S_, "saturation_liquid", tag_next_, tags_[1].second);
  // aliasVector(*S_, "porosity", tag_next_, tags_[1].second);
  // aliasVector(*S_, "molar_density_liquid", tag_next_, tags_[1].second);
  // aliasVector(*S_, "surface-water_flux", tag_next_, tags_[1].second);
  // aliasVector(*S_, "surface-ponded_depth", tag_next_, tags_[1].second);
  // aliasVector(*S_, "surface-porosity", tag_next_, tags_[1].second);
  // aliasVector(*S_, "surface-molar_density_liquid", tag_next_, tags_[1].second);

  sub_pks_[1]->Setup();
  // END HACK
}

void
MPCSubcycled::Initialize()

{
  int i = 0;
  for (const auto& tag : tags_) {
    if (subcycling_[i])
      S_->GetRecordW("dt", tag.second, name()).set_initialized();
    ++i;
  }
  MPC<PK>::Initialize();
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCSubcycled::get_dt()
{
  double dt = std::numeric_limits<double>::max();
  int i = 0;
  for (auto& pk : sub_pks_) {
    dts_[i] = pk->get_dt();
    if (!subcycling_[i]) dt = std::min(dt, dts_[i]);
    ++i;
  }

  for (auto& dt_local : dts_) dt_local = std::min(dt_local, dt);
  dt_ = dt;
  return dt;
}


// -----------------------------------------------------------------------------
// Set standard dt
// -----------------------------------------------------------------------------
void MPCSubcycled::set_dt(double dt)
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
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool MPCSubcycled::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  AMANZI_ASSERT(std::abs(t_new - t_old - dt_) < 1.e-4);

  int i = 0;
  for (auto& pk : sub_pks_) {
    if (subcycling_[i]) {
      // advance the subcycled, subcycling if needed
      bool done = false;
      double t_inner = t_old;
      double dt_inner = dts_[i];

      Tag tag_subcycle_current = tags_[i].first;
      Tag tag_subcycle_next = tags_[i].second;

      S_->set_time(tag_subcycle_current, t_old);
      while (!done) {
        dt_inner = std::min(dt_inner, t_new - t_inner);
        S_->Assign("dt", tag_subcycle_next, name(), dt_inner);
        S_->set_time(tag_subcycle_next, t_inner + dt_inner);
        bool fail_inner = sub_pks_[i]->AdvanceStep(t_inner, t_inner+dt_inner, false);

        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  step failed? " << fail_inner << std::endl;
        bool valid_inner = sub_pks_[i]->ValidStep();
        if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
          *vo_->os() << "  step valid? " << valid_inner << std::endl;
        }

        if (fail_inner || !valid_inner) {
          sub_pks_[i]->FailStep(t_old, t_new, tag_subcycle_next);

          dt_inner = sub_pks_[i]->get_dt();
          S_->set_time(tag_subcycle_next, S_->get_time(tag_subcycle_current));

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

        if (dt_inner < min_dt_) {
          Errors::Message msg;
          msg << "SubPK " << pk->name() << " crashing timestep in subcycling: dt = " << dt_inner;
          Exceptions::amanzi_throw(msg);
        }
      }

    } else {
      // advance the standard PK using the full step size
      bool fail = pk->AdvanceStep(t_old, t_new, reinit);
      if (fail) return fail;
    }
    ++i;
  }

  return false;
}


// void
// MPCSubcycled::CommitStep(double t_old, double t_new, const Tag& tag)
// {
//   // do not commitstep on subcycled -- this has already been done
//   int i = 0;
//   for (auto& pk : sub_pks_) {
//     if (!subcycling_[i]) pk->CommitStep(t_old, t_new, tag);
//     ++i;
//   }
// }



}  // namespace Amanzi

