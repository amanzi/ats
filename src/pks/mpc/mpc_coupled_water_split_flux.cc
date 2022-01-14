/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "EvaluatorPrimary.hh"
#include "mpc_surface_subsurface_helpers.hh"

#include "mpc_coupled_water_split_flux.hh"

#include "PK_Physical.hh"

namespace Amanzi {

MPCCoupledWaterSplitFlux::MPCCoupledWaterSplitFlux(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPC<PK>(FElist, plist, S, solution),
      eval_pvfe_(Teuchos::null)
{
  // collect keys and names
  std::string domain = plist_->get<std::string>("domain name");
  std::string domain_star = plist_->get<std::string>("star domain name", domain+"_star");
  primary_variable_ = Keys::readKey(*plist_, domain, "primary variable");
  primary_variable_star_ = Keys::readKey(*plist_, domain_star, "primary variable star", Keys::getVarName(primary_variable_));
  lateral_flow_source_ = Keys::readKey(*plist_, domain, "lateral flow source", "lateral_flow_source");
  conserved_variable_star_ = Keys::readKey(*plist_, domain_star, "conserved quantity star", "water_content");
  cv_key_ = Keys::readKey(*plist_, domain, "cell volume", "cell_volume");
  
  // set up for a primary variable field evaluator for the flux
  auto& sublist = S->FEList().sublist(lateral_flow_source_);
  sublist.set("field evaluator type", "primary variable");

  // init sub-pks
  init_(S);
};


// -- initialize in reverse order
void MPCCoupledWaterSplitFlux::Initialize(const Teuchos::Ptr<State>& S)
{
  sub_pks_[1]->Initialize(S);
  CopyPrimaryToStar(S, S);
  S->GetField(primary_variable_star_, "star")->set_initialized();
  sub_pks_[0]->Initialize(S);
}


void MPCCoupledWaterSplitFlux::Setup(const Teuchos::Ptr<State>& S) {
  MPC<PK>::Setup(S);

  S->RequireEvaluator(lateral_flow_source_);
}

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCCoupledWaterSplitFlux::get_dt() {
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
void MPCCoupledWaterSplitFlux::set_dt( double dt) {
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }
};

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool MPCCoupledWaterSplitFlux::AdvanceStep(double t_old, double t_new, bool reinit) {
  // Advance the star system 
  bool fail = false;
  AMANZI_ASSERT(sub_pks_.size() == 2);
  fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // Copy star's new value into primary's old value
  CopyStarToPrimary(t_new - t_old);

  // Now advance the primary
  fail = sub_pks_[1]->AdvanceStep(t_old, t_new, reinit);
  return fail;
};


void MPCCoupledWaterSplitFlux::CommitStep(double t_old, double t_new,
        const Teuchos::RCP<State>& S) {
  // commit before copy to ensure record for extrapolation in star system uses
  // its own solutions
  MPC<PK>::CommitStep(t_old, t_new, S);

  // Copy the primary into the star to advance
  CopyPrimaryToStar(S.ptr(), S.ptr());
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyPrimaryToStar(const Teuchos::Ptr<const State>& S,
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
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyStarToPrimary(double dt) {
  // make sure we have the evaluator at the new state timestep
  if (eval_pvfe_ == Teuchos::null) {
    Teuchos::RCP<Evaluator> fe = S_next_->GetEvaluator(lateral_flow_source_);
    eval_pvfe_ = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(fe);
    AMANZI_ASSERT(eval_pvfe_ != Teuchos::null);
  }

  auto star_pk = Teuchos::rcp_dynamic_cast<PK_Physical>(sub_pks_[0]);
  AMANZI_ASSERT(star_pk.get());

  // these updates should do nothing, but you never know
  S_inter_->GetEvaluator(conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetEvaluator(conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);

  // grab the data, difference
  star_pk->debugger()->WriteVector("WC0", S_inter_->GetPtr<CompositeVector>(conserved_variable_star_).ptr());
  star_pk->debugger()->WriteVector("WC1", S_next_->GetPtr<CompositeVector>(conserved_variable_star_).ptr());
  auto& q_div = *S_next_->GetW<CompositeVector>(lateral_flow_source_, S_next_->GetField(lateral_flow_source_).owner())
                ->ViewComponent("cell",false);
  q_div.Update(1.0/dt,
               *S_next_->Get<CompositeVector>(conserved_variable_star_).ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->Get<CompositeVector>(conserved_variable_star_).ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_next_->Get<CompositeVector>(cv_key_).ViewComponent("cell",false), q_div, 0.);
  star_pk->debugger()->WriteVector("qdiv", S_next_->GetPtr<CompositeVector>(lateral_flow_source_).ptr());

  // mark the source evaluator as changed to ensure the total source gets updated.
  eval_pvfe_->SetChanged(S_next_.ptr());
}

} // namespace
