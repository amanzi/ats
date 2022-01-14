/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "EvaluatorPrimary.hh"
#include "mpc_surface_subsurface_helpers.hh"

#include "mpc_permafrost_split_flux.hh"

#include "PK_Physical.hh"

namespace Amanzi {

MPCPermafrostSplitFlux::MPCPermafrostSplitFlux(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPC<PK>(FElist, plist, S, solution),
      p_eval_pvfe_(Teuchos::null),
      T_eval_pvfe_(Teuchos::null)
{
  // collect keys and names
  std::string domain = plist_->get<std::string>("domain name");
  std::string domain_star = plist_->get<std::string>("star domain name", domain+"_star");
  p_primary_variable_ = Keys::readKey(*plist_, domain, "pressure primary variable", "pressure");
  T_primary_variable_ = Keys::readKey(*plist_, domain, "temperature primary variable", "temperature");

  p_primary_variable_star_ = Keys::readKey(*plist_, domain_star, "pressure primary variable star", Keys::getVarName(p_primary_variable_));
  T_primary_variable_star_ = Keys::readKey(*plist_, domain_star, "temperature primary variable star", Keys::getVarName(T_primary_variable_));

  p_lateral_flow_source_ = Keys::readKey(*plist_, domain, "mass lateral flow source", "mass_lateral_flow_source");
  T_lateral_flow_source_ = Keys::readKey(*plist_, domain, "energy lateral flow source", "energy_lateral_flow_source");

  p_conserved_variable_star_ = Keys::readKey(*plist_, domain_star, "mass conserved quantity star", "water_content");
  T_conserved_variable_star_ = Keys::readKey(*plist_, domain_star, "energy conserved quantity star", "energy");

  cv_key_ = Keys::readKey(*plist_, domain, "cell volume", "cell_volume");
  
  // set up for a primary variable field evaluator for the flux
  auto& p_sublist = S->FEList().sublist(p_lateral_flow_source_);
  p_sublist.set("field evaluator type", "primary variable");
  auto& T_sublist = S->FEList().sublist(T_lateral_flow_source_);
  T_sublist.set("field evaluator type", "primary variable");

  // init sub-pks
  init_(S);
};


void MPCPermafrostSplitFlux::Initialize(const Teuchos::Ptr<State>& S)
{
  sub_pks_[1]->Initialize(S);
  CopyPrimaryToStar(S, S);
  S->GetField(p_primary_variable_star_, S->GetField(p_primary_variable_star_)->owner())->set_initialized();
  S->GetField(T_primary_variable_star_, S->GetField(T_primary_variable_star_)->owner())->set_initialized();
  sub_pks_[0]->Initialize(S);
}


void MPCPermafrostSplitFlux::Setup(const Teuchos::Ptr<State>& S)
{
  MPC<PK>::Setup(S);

  S->Require<CompositeVector,CompositeVectorSpace>(p_lateral_flow_source_, Tags::NEXT,  p_lateral_flow_source_)
      ->SetMesh(S->GetMesh(Keys::getDomain(p_lateral_flow_source_)))
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(p_lateral_flow_source_);

  S->Require<CompositeVector,CompositeVectorSpace>(T_lateral_flow_source_, Tags::NEXT,  T_lateral_flow_source_)
      ->SetMesh(S->GetMesh(Keys::getDomain(T_lateral_flow_source_)))
      ->SetComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(T_lateral_flow_source_);
}

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCPermafrostSplitFlux::get_dt()
{
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
void MPCPermafrostSplitFlux::set_dt( double dt)
{
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }
};

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool MPCPermafrostSplitFlux::AdvanceStep(double t_old, double t_new, bool reinit)
{
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


void MPCPermafrostSplitFlux::CommitStep(double t_old, double t_new,
        const Teuchos::RCP<State>& S)
{
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
MPCPermafrostSplitFlux::CopyPrimaryToStar(const Teuchos::Ptr<const State>& S,
                                    const Teuchos::Ptr<State>& S_star)
{
  // copy p primary variable
  auto& p_star = *S_star->GetW<CompositeVector>(p_primary_variable_star_, S_star->GetField(p_primary_variable_star_).owner())
                  ->ViewComponent("cell",false);
  const auto& p = *S->Get<CompositeVector>(p_primary_variable_).ViewComponent("cell",false);
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p[0][c] <= 101325.0) {
      p_star[0][c] = 101325.;
    } else {
      p_star[0][c] = p[0][c];
    }
  }

  auto peval = S_star->GetEvaluator(p_primary_variable_star_);
  auto peval_pvfe = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(peval);
  peval_pvfe->SetChanged(S_star.ptr());

  // copy T primary variable
  auto& T_star = *S_star->GetW<CompositeVector>(T_primary_variable_star_, S_star->GetField(T_primary_variable_star_).owner())
                  ->ViewComponent("cell",false);
  const auto& T = *S->Get<CompositeVector>(T_primary_variable_).ViewComponent("cell",false);
  T_star = T;

  auto Teval = S_star->GetEvaluator(T_primary_variable_star_);
  auto Teval_pvfe = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(Teval);
  Teval_pvfe->SetChanged(S_star.ptr());

}

// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary(double dt)
{
  // make sure we have the evaluator at the new state timestep
  if (p_eval_pvfe_ == Teuchos::null) {
    Teuchos::RCP<Evaluator> fe = S_next_->GetEvaluator(p_lateral_flow_source_);
    p_eval_pvfe_ = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(fe);
    AMANZI_ASSERT(p_eval_pvfe_ != Teuchos::null);
  }
  
  if (T_eval_pvfe_ == Teuchos::null) {
    Teuchos::RCP<Evaluator> fe = S_next_->GetEvaluator(T_lateral_flow_source_);
    T_eval_pvfe_ = Teuchos::rcp_dynamic_cast<EvaluatorPrimary>(fe);
    AMANZI_ASSERT(T_eval_pvfe_ != Teuchos::null);
  }

  // these updates should do nothing, but you never know
  S_inter_->GetEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);

  // grab the data, difference
  auto& q_div = *S_next_->GetW<CompositeVector>(p_lateral_flow_source_, S_next_->GetField(p_lateral_flow_source_).owner())
                ->ViewComponent("cell",false);
  q_div.Update(1.0/dt,
               *S_next_->Get<CompositeVector>(p_conserved_variable_star_).ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->Get<CompositeVector>(p_conserved_variable_star_).ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_next_->Get<CompositeVector>(cv_key_).ViewComponent("cell",false), q_div, 0.);

  // mark the source evaluator as changed to ensure the total source gets updated.
  p_eval_pvfe_->SetChanged(S_next_.ptr());

  // grab the data, difference
  auto& qE_div = *S_next_->GetW<CompositeVector>(T_lateral_flow_source_, S_next_->GetField(T_lateral_flow_source_).owner())
                ->ViewComponent("cell",false);
  qE_div.Update(1.0/dt,
               *S_next_->Get<CompositeVector>(T_conserved_variable_star_).ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->Get<CompositeVector>(T_conserved_variable_star_).ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_next_->Get<CompositeVector>(cv_key_).ViewComponent("cell",false), qE_div, 0.);

  // mark the source evaluator as changed to ensure the total source gets updated.
  T_eval_pvfe_->SetChanged(S_next_.ptr());
}

} // namespace
