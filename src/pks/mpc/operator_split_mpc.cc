/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived OperatorSplitMPC class.  Provides only the advance()
method missing from MPC.hh.  In operatorsplit coupling, we simply loop over the
sub-PKs, calling their advance() methods and returning failure if any fail.

NOTE this is not general, and shouldn't be used arbitrarily.  Need to
reimplement for general case.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include "EvaluatorPrimary.hh"
#include "mpc_surface_subsurface_helpers.hh"

#include "operator_split_mpc.hh"


namespace Amanzi {

OperatorSplitMPC::OperatorSplitMPC(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPC<PK>(FElist, plist, S, solution) {

  std::string domain = plist_->get<std::string>("domain name");
  primary_variable_ = Keys::readKey(*plist_, domain, "primary variable");
  primary_variable_star_ = Keys::getKey(domain+"_star", Keys::getVarName(primary_variable_));
  init_(S);
};


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double OperatorSplitMPC::get_dt() {
  double dt = 1.0e99;
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt, (*pk)->get_dt());
  }
  return dt;
};

// -----------------------------------------------------------------------------
// Set timestep for sub PKs 
// -----------------------------------------------------------------------------
void OperatorSplitMPC::set_dt( double dt) {
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }

};

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool OperatorSplitMPC::AdvanceStep(double t_old, double t_new, bool reinit) {

  // Advance the star system 
  bool fail = false;
  AMANZI_ASSERT(sub_pks_.size() == 2);
  fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // Copy star's new value into primary's old value
  CopyStarToPrimary(S_next_.ptr(), S_inter_.ptr());
  CopyStarToPrimary(S_next_.ptr(), S_next_.ptr());

  // BEGIN THE NON-GENERIC PART TO BE REMOVED
  // also copy and mark the subsurface system
  CopySurfaceToSubsurface(*S_inter_->GetPtr<CompositeVector>(primary_variable_),
                          S_inter_->GetW<CompositeVector>("pressure",S_inter_->GetField("pressure").owner()).ptr());
  auto eval = S_inter_->GetEvaluator("pressure");
  auto eval_pvfe = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(eval);
  eval_pvfe->SetChanged(S_inter_.ptr());

  CopySurfaceToSubsurface(*S_next_->GetPtr<CompositeVector>(primary_variable_),
                          S_next_->GetW<CompositeVector>("pressure",S_next_->GetField("pressure").owner()).ptr());
  auto eval2 = S_next_->GetEvaluator("pressure");
  auto eval_pvfe2 = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(eval2);
  eval_pvfe2->SetChanged(S_next_.ptr());
  // END THE NON-GENERIC PART TO BE REMOVED


  // Now advance the primary
  fail = sub_pks_[1]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // Now copy back to star.
  CopyPrimaryToStar(S_next_.ptr(), S_next_.ptr());
  return fail;
};


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
OperatorSplitMPC::CopyPrimaryToStar(const Teuchos::Ptr<const State>& S,
                                    const Teuchos::Ptr<State>& S_star) {
  auto& pv_star = *S_star->GetW<CompositeVector>(primary_variable_star_, S_star->GetField(primary_variable_star_).owner())
                  ->ViewComponent("cell",false);
  auto& pv = *S->GetPtr<CompositeVector>(primary_variable_)
             ->ViewComponent("cell",false);
  for (int c=0; c!=pv_star.MyLength(); ++c) {
    if (pv[0][c] <= 101325.0) {
      pv_star[0][c] = 101325.;
    } else {
      pv_star[0][c] = pv[0][c];
    }
  }

  auto eval = S_star->GetEvaluator(primary_variable_star_);
  auto eval_pvfe = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(eval);
  eval_pvfe->SetChanged(S_star.ptr());
}

// -----------------------------------------------------------------------------
// Copy the star variable to the primary system
// -----------------------------------------------------------------------------
void
OperatorSplitMPC::CopyStarToPrimary(const Teuchos::Ptr<const State>& S_star,
                                    const Teuchos::Ptr<State>& S) {
  auto& pv_star = *S_star->GetPtr<CompositeVector>(primary_variable_star_)
                  ->ViewComponent("cell",false);
  auto& pv = *S->GetW<CompositeVector>(primary_variable_, S->GetField(primary_variable_).owner())
                  ->ViewComponent("cell",false);
  for (int c=0; c!=pv_star.MyLength(); ++c) {
    if (pv_star[0][c] > 101325.0) {
      pv[0][c] = pv_star[0][c];
    }
  }

  auto eval = S->GetEvaluator(primary_variable_);
  auto eval_pvfe = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(eval);
  eval_pvfe->SetChanged(S.ptr());
}


} // namespace
