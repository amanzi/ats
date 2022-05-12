/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "mpc_coupled_water_split_flux.hh"

#include "mpc_surface_subsurface_helpers.hh"
#include "pk_helpers.hh"
#include "PK_Physical.hh"

namespace Amanzi {

MPCCoupledWaterSplitFlux::MPCCoupledWaterSplitFlux(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPC<PK>(FElist, plist, S, solution)
{
  // collect keys and names
  std::string domain = Keys::readDomain(*plist_);
  std::string domain_star = Keys::readDomain(*plist_, "star", domain+"_star");
  primary_variable_ = Keys::readKey(*plist_, domain, "primary variable");
  primary_variable_star_ = Keys::readKey(*plist_, domain_star, "primary variable star", Keys::getVarName(primary_variable_));
  lateral_flow_source_ = Keys::readKey(*plist_, domain, "lateral flow source", "lateral_flow_source");
  conserved_variable_star_ = Keys::readKey(*plist_, domain_star, "conserved quantity star", "water_content");
  cv_key_ = Keys::readKey(*plist_, domain, "cell volume", "cell_volume");

  // init sub-pks
  init_();
  AMANZI_ASSERT(sub_pks_.size() == 2);
};


// -- initialize in reverse order
void MPCCoupledWaterSplitFlux::Initialize()
{
  sub_pks_[1]->Initialize();
  CopyPrimaryToStar(tag_next_, tag_next_);
  S_->GetRecordW(primary_variable_star_, tag_next_, primary_variable_star_).set_initialized();
  sub_pks_[0]->Initialize();
}


void MPCCoupledWaterSplitFlux::Setup()
{
  MPC<PK>::Setup();
  S_->Require<CompositeVector,CompositeVectorSpace>(lateral_flow_source_, tag_next_);
  RequireEvaluatorPrimary(lateral_flow_source_, tag_next_, *S_);
}

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCCoupledWaterSplitFlux::get_dt()
{
  double dt = std::numeric_limits<double>::max();
  for (auto& pk : sub_pks_) dt = std::min(dt, pk->get_dt());
  return dt;
};

// -----------------------------------------------------------------------------
// Set timestep for sub PKs
// -----------------------------------------------------------------------------
void MPCCoupledWaterSplitFlux::set_dt( double dt)
{
  for (auto& pk : sub_pks_) pk->set_dt(dt);
};

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool MPCCoupledWaterSplitFlux::AdvanceStep(double t_old, double t_new, bool reinit)
{
  // Advance the star system
  bool fail = false;
  fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // Copy star's new value into primary's old value
  CopyStarToPrimary(tag_current_, tag_next_, tag_current_);

  // Now advance the primary
  fail = sub_pks_[1]->AdvanceStep(t_old, t_new, reinit);
  return fail;
};


void MPCCoupledWaterSplitFlux::CommitStep(double t_old, double t_new,
        const Tag& tag)
{
  // commit before copy to ensure record for extrapolation in star system uses
  // its own solutions
  MPC<PK>::CommitStep(t_old, t_new, tag);

  // Copy the primary into the star to advance
  CopyPrimaryToStar(tag, tag);
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyPrimaryToStar(const Tag& tag_primary, const Tag& tag_star)
{
  auto& pv_star = *S_->GetW<CompositeVector>(primary_variable_star_, tag_star, primary_variable_star_)
                  .ViewComponent("cell",false);
  auto& pv = *S_->Get<CompositeVector>(primary_variable_, tag_primary)
                  .ViewComponent("cell",false);
  for (int c=0; c!=pv_star.MyLength(); ++c) {
    if (pv[0][c] <= 101325.0) {
      pv_star[0][c] = 101325.;
    } else {
      pv_star[0][c] = pv[0][c];
    }
  }

  ChangedEvaluatorPrimary(primary_variable_star_, tag_star, *S_);
}

// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyStarToPrimary(const Tag& tag_star_current, const Tag& tag_star_next,
        const Tag& tag_primary)
{
  auto star_pk = Teuchos::rcp_dynamic_cast<PK_Physical>(sub_pks_[0]);
  AMANZI_ASSERT(star_pk.get());

  // these updates should do nothing, but you never know
  S_->GetEvaluator(conserved_variable_star_, tag_star_current).Update(*S_, name_);
  S_->GetEvaluator(conserved_variable_star_, tag_star_next).Update(*S_, name_);

  // grab the data, difference
  star_pk->debugger()->WriteVector("WC0", S_->GetPtr<CompositeVector>(conserved_variable_star_, tag_star_current).ptr());
  star_pk->debugger()->WriteVector("WC1", S_->GetPtr<CompositeVector>(conserved_variable_star_, tag_star_next).ptr());
  double dt = S_->get_time(tag_star_next) - S_->get_time(tag_star_current);
  auto& q_div = *S_->GetW<CompositeVector>(lateral_flow_source_, tag_primary, lateral_flow_source_)
                .ViewComponent("cell",false);
  q_div.Update(1.0/dt,
               *S_->Get<CompositeVector>(conserved_variable_star_, tag_star_next).ViewComponent("cell",false),
               -1.0/dt,
               *S_->Get<CompositeVector>(conserved_variable_star_, tag_star_current).ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_->Get<CompositeVector>(cv_key_, tag_primary).ViewComponent("cell",false), q_div, 0.);
  star_pk->debugger()->WriteVector("qdiv", S_->GetPtr<CompositeVector>(lateral_flow_source_, tag_primary).ptr());

  // mark the source evaluator as changed to ensure the total source gets updated.
  ChangedEvaluatorPrimary(lateral_flow_source_, tag_primary, *S_);
}

} // namespace
