/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.


*/

/*

*/

#include "mpc_urban_hydro_subcycled.hh"
#include "pk_helpers.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCUrbanHydroSub::MPCUrbanHydroSub(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln) :
  MPCSubcycled(pk_tree, global_list, S, soln)
{
  pipe_domain_ = "network";
  master_domain_ = "surface";

  pipe_flow_src_key_ = Keys::readKey(*plist_, pipe_domain_, "pipe source", "source_drain"); ;
  master_src_key_ = Keys::readKey(*plist_, master_domain_, "master", "source_drain"); ;

  master_id_ = 0;
  subcycled_id_ = 1;

    
}


// -----------------------------------------------------------------------------
// Advance the ith sub-PK individually, returning a failure
// -----------------------------------------------------------------------------
// bool
// MPCUrbanHydroSub::AdvanceStep_i_(std::size_t i, double t_old, double t_new, bool reinit)
// {
//   bool fail = false;
//   if (subcycling_[i]) {
//     // advance the subcycled, subcycling if needed
//     bool done = false;
//     double t_inner = t_old;
//     double dt_inner = dts_[i];

//     Tag tag_subcycle_current = tags_[i].first;
//     Tag tag_subcycle_next = tags_[i].second;

//     S_->set_time(tag_subcycle_current, t_old);
//     while (!done) {
//       dt_inner = std::min(dt_inner, t_new - t_inner);
//       S_->Assign("dt", tag_subcycle_current, name(), dt_inner);
//       S_->set_time(tag_subcycle_next, t_inner + dt_inner);
//       bool fail_inner = sub_pks_[i]->AdvanceStep(t_inner, t_inner+dt_inner, false);

//       if (vo_->os_OK(Teuchos::VERB_EXTREME))
//         *vo_->os() << "  step failed? " << fail_inner << std::endl;
//       bool valid_inner = sub_pks_[i]->ValidStep();
//       if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
//         *vo_->os() << "  step valid? " << valid_inner << std::endl;
//       }

//       if (fail_inner || !valid_inner) {
//         sub_pks_[i]->FailStep(t_old, t_new, tag_subcycle_next);

//         dt_inner = sub_pks_[i]->get_dt();
//         S_->set_time(tag_subcycle_next, S_->get_time(tag_subcycle_current));

//         if (vo_->os_OK(Teuchos::VERB_EXTREME))
//           *vo_->os() << "  failed, new timestep is " << dt_inner << std::endl;

//       } else {
//         sub_pks_[i]->CommitStep(t_inner, t_inner + dt_inner, tag_subcycle_next);
//         t_inner += dt_inner;
//         if (std::abs(t_new - t_inner) < 1.e-10) done = true;

//         S_->set_time(tag_subcycle_current, S_->get_time(tag_subcycle_next));
//         S_->advance_cycle(tag_subcycle_next);

//         dt_inner = sub_pks_[i]->get_dt();
//         if (vo_->os_OK(Teuchos::VERB_EXTREME))
//           *vo_->os() << "  success, new timestep is " << dt_inner << std::endl;
//       }
//     }

//   } else {
//     // advance the standard PK using the full step size
//     fail = sub_pks_[i]->AdvanceStep(t_old, t_new, reinit);
//   }
//   return fail;
// }



// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool MPCUrbanHydroSub::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;

  bool done = false;
  double t_inner = t_old;
  double dt_inner = dts_[subcycled_id_];

  Tag tag_subcycle_current = tags_[subcycled_id_].first;
  Tag tag_subcycle_next = tags_[subcycled_id_].second;

  S_->set_time(tag_subcycle_current, t_old);
  while (!done) {
    dt_inner = std::min(dt_inner, t_new - t_inner);
    S_->Assign("dt", tag_subcycle_current, name(), dt_inner);
    S_->set_time(tag_subcycle_next, t_inner + dt_inner);
    bool fail_inner = sub_pks_[subcycled_id_]->AdvanceStep(t_inner, t_inner+dt_inner, false);

    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "  step failed? " << fail_inner << std::endl;
    bool valid_inner = sub_pks_[subcycled_id_]->ValidStep();
    if (vo_->os_OK(Teuchos::VERB_EXTREME)) {
      *vo_->os() << "  step valid? " << valid_inner << std::endl;
    }

    if (fail_inner || !valid_inner) {
      sub_pks_[subcycled_id_]->FailStep(t_old, t_new, tag_subcycle_next);
      
      dt_inner = sub_pks_[subcycled_id_]->get_dt();
      S_->set_time(tag_subcycle_next, S_->get_time(tag_subcycle_current));
      
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  failed, new timestep is " << dt_inner << std::endl;

    } else {
      sub_pks_[subcycled_id_]->CommitStep(t_inner, t_inner + dt_inner, tag_subcycle_next);
      t_inner += dt_inner;
      if (std::abs(t_new - t_inner) < 1.e-10) done = true;
      
      S_->set_time(tag_subcycle_current, S_->get_time(tag_subcycle_next));
      S_->advance_cycle(tag_subcycle_next);
      
      dt_inner = sub_pks_[subcycled_id_]->get_dt();
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  success, new timestep is " << dt_inner << std::endl;
    }
  }

  
  int i = 0;
  for (auto& pk : sub_pks_) {
    fail = AdvanceStep_i_(i, t_old, t_new, reinit);
    if (fail) return fail;
    ++i;
  }

  return fail;
}


// void
// MPCUrbanHydroSub::CommitStep(double t_old, double t_new, const Tag& tag)
// {
//   MPC<PK>::CommitStep(t_old, t_new, tag);

//   if (S_->get_cycle() < 0 && tag == Tags::NEXT) {
//     // initial commit, also do the substep commits
//     int i = 0;
//     for (auto& pk : sub_pks_) {
//       if (subcycling_[i]) {
//         pk->CommitStep(t_old, t_new, tags_[i].second);
//       }
//       ++i;
//     }
//   }
// }



}  // namespace Amanzi

