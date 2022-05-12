/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "mpc_permafrost_split_flux.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "pk_helpers.hh"

namespace Amanzi {

MPCPermafrostSplitFlux::MPCPermafrostSplitFlux(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPCSubcycled(FElist, plist, S, solution)
{
  // collect domain names
  domain_set_ = Keys::readDomain(*plist_); // e.g. surface or surface_column:*
  domain_star_ = Keys::readDomain(*plist_, "star"); // e.g. surface_star

  // determine whether we are coupling subdomains or coupling 3D domains
  is_domain_set_ = S_->HasDomainSet(domain_set_);
  if (is_domain_set_) domain_ = Keys::getDomainInSet(domain_set_, "*");
  else domain_ = domain_set_;

  domain_sub_ = Keys::readDomainHint(*plist_, domain_set_, "surface", "subsurface");
  domain_snow_ = Keys::readDomainHint(*plist_, domain_set_, "surface", "snow");

  // determine the coupling strategy: "pressure" passes the pressure field,
  // "flux" the flux field, while "hybrid" passes one or the other depending
  // upon conditions.  "hybrid" is the most robust.
  coupling_ = plist_->get<std::string>("coupling type", "hybrid");
  if (coupling_ != "pressure" && coupling_ != "flux" && coupling_ != "hybrid") {
    Errors::Message msg("WeakMPCSemiCoupled: \"coupling type\" must be one of \"pressure\", \"flux\", or \"hybrid\".");
    Exceptions::amanzi_throw(msg);
  }

  // collect keys and names
  // -- primary variable for the main surface domain
  p_primary_variable_ = Keys::readKey(*plist_, domain_, "pressure primary variable", "pressure");
  T_primary_variable_ = Keys::readKey(*plist_, domain_, "temperature primary variable", "temperature");
  p_primary_variable_suffix_ = Keys::getVarName(p_primary_variable_);
  T_primary_variable_suffix_ = Keys::getVarName(T_primary_variable_);

  // -- primary variable for the main subsurface domain
  p_sub_primary_variable_ = Keys::readKey(*plist_, domain_sub_, "subsurface pressure primary variable", "pressure");
  T_sub_primary_variable_ = Keys::readKey(*plist_, domain_sub_, "subsurface temperature primary variable", "temperature");
  p_sub_primary_variable_suffix_ = Keys::getVarName(p_sub_primary_variable_);
  T_sub_primary_variable_suffix_ = Keys::getVarName(T_sub_primary_variable_);

  // -- need to save the updated conserved quantity too
  p_conserved_variable_ = Keys::readKey(*plist_, domain_, "water conserved quantity", "water_content");
  T_conserved_variable_ = Keys::readKey(*plist_, domain_, "energy conserved quantity", "energy");
  p_conserved_variable_star_ = Keys::readKey(*plist_, domain_star_, "water conserved quantity star", "water_content");
  T_conserved_variable_star_ = Keys::readKey(*plist_, domain_star_, "energy conserved quantity star", "energy");

  // -- primary variable for the star domain
  p_primary_variable_star_ = Keys::readKey(*plist_, domain_star_,
          "pressure primary variable star", Keys::getVarName(p_primary_variable_));
  T_primary_variable_star_ = Keys::readKey(*plist_, domain_star_,
          "temperature primary variable star", Keys::getVarName(T_primary_variable_));

  // -- flux variables for coupling
  if (coupling_ != "pressure") {
    p_lateral_flow_source_ = Keys::readKey(*plist_, domain_, "water lateral flow source", "water_lateral_flow_source");
    p_lateral_flow_source_suffix_ = Keys::getVarName(p_lateral_flow_source_);
    T_lateral_flow_source_ = Keys::readKey(*plist_, domain_, "energy lateral flow source", "energy_lateral_flow_source");
    T_lateral_flow_source_suffix_ = Keys::getVarName(T_lateral_flow_source_);

    cv_key_ = Keys::readKey(*plist_, domain_star_, "cell volume", "cell_volume");
  }
};


void MPCPermafrostSplitFlux::Setup()
{
  MPCSubcycled::Setup();

  if (coupling_ != "pressure") {
    if (is_domain_set_) {
      auto domain_set = S_->GetDomainSet(domain_set_);
      for (const auto& domain : *domain_set) {
        auto p_key = Keys::getKey(domain, p_lateral_flow_source_suffix_);
        S_->Require<CompositeVector,CompositeVectorSpace>(p_key, tag_next_, p_key)
          .SetMesh(S_->GetMesh(domain))
          ->SetComponent("cell", AmanziMesh::CELL, 1);
        RequireEvaluatorPrimary(p_key, tag_next_, *S_);

        auto T_key = Keys::getKey(domain, T_lateral_flow_source_suffix_);
        S_->Require<CompositeVector,CompositeVectorSpace>(T_key, tag_next_, T_key)
          .SetMesh(S_->GetMesh(domain))
          ->SetComponent("cell", AmanziMesh::CELL, 1);
        RequireEvaluatorPrimary(T_key, tag_next_, *S_);
      }
    } else {
      S_->Require<CompositeVector,CompositeVectorSpace>(p_lateral_flow_source_, tag_next_,  p_lateral_flow_source_)
        .SetMesh(S_->GetMesh(Keys::getDomain(p_lateral_flow_source_)))
        ->SetComponent("cell", AmanziMesh::CELL, 1);
      RequireEvaluatorPrimary(p_lateral_flow_source_, tag_next_, *S_);

      S_->Require<CompositeVector,CompositeVectorSpace>(T_lateral_flow_source_, tag_next_,  T_lateral_flow_source_)
        .SetMesh(S_->GetMesh(Keys::getDomain(T_lateral_flow_source_)))
        ->SetComponent("cell", AmanziMesh::CELL, 1);
      RequireEvaluatorPrimary(T_lateral_flow_source_, tag_next_, *S_);
    }

    // also need conserved quantities at old and new times
    S_->Require<CompositeVector,CompositeVectorSpace>(p_conserved_variable_star_, tag_next_)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->Require<CompositeVector,CompositeVectorSpace>(p_conserved_variable_star_, tag_current_)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(p_conserved_variable_star_, tag_next_);
    //S_->RequireEvaluator(p_conserved_variable_star_, tag_current_);

    S_->Require<CompositeVector,CompositeVectorSpace>(T_conserved_variable_star_, tag_next_)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->Require<CompositeVector,CompositeVectorSpace>(T_conserved_variable_star_, tag_current_)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S_->RequireEvaluator(T_conserved_variable_star_, tag_next_);
    //S_->RequireEvaluator(T_conserved_variable_star_, tag_current_);
  }
}



void MPCPermafrostSplitFlux::Initialize()
{
  sub_pks_[1]->Initialize();
  CopyPrimaryToStar_(tag_next_, tag_next_);
  sub_pks_[0]->Initialize();

  // FIXME -- this order is logically wrong but currently necessary.  The
  // intention is the following initialization process:
  //
  // 1. CopyPrimaryToStar() sets surface_star cell values of p & T
  // 2. sub_pk->Initialize() sets face values from cell values.
  //
  // Logically, one would think that setting these as initialized before
  // calling sub_pk->Initialize() is right then -- since cell values are
  // already set, and setting face values is "non-initializing", this seems
  // right.  But for some reason, that bypasses the call to
  // SetFaceValuesFromCells() or whatever.  --etc
  auto p_owner = S_->GetRecord(p_primary_variable_star_, tag_next_).owner();
  S_->GetRecordW(p_primary_variable_star_, tag_next_, p_owner).set_initialized();
  auto T_owner = S_->GetRecord(T_primary_variable_star_, tag_next_).owner();
  S_->GetRecordW(T_primary_variable_star_, tag_next_, T_owner).set_initialized();

  auto& pstar = *S_->Get<CompositeVector>(p_primary_variable_star_, tag_next_)
    .ViewComponent("cell", false);

  // set the fluxes as initialized -- they will get set by later calls to
  // CopyStarToPrimary, so no need for values, but do need to toggle the flag.
  if (coupling_ != "pressure") {
    if (is_domain_set_) {
      auto domain_set = S_->GetDomainSet(domain_set_);
      for (const auto& domain : *domain_set) {
        auto pkey = Keys::getKey(domain, p_lateral_flow_source_suffix_);
        auto p_owner = S_->GetRecord(pkey, tag_next_).owner();
        S_->GetRecordW(pkey, tag_next_, p_owner).set_initialized();

        auto Tkey = Keys::getKey(domain, T_lateral_flow_source_suffix_);
        auto T_owner = S_->GetRecord(Tkey, tag_next_).owner();
        S_->GetRecordW(Tkey, tag_next_, T_owner).set_initialized();
      }
    } else {
      S_->GetRecordW(p_lateral_flow_source_, tag_next_, p_lateral_flow_source_).set_initialized();
      S_->GetRecordW(T_lateral_flow_source_, tag_next_, T_lateral_flow_source_).set_initialized();
    }
  }
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool MPCPermafrostSplitFlux::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  AMANZI_ASSERT(std::abs(t_new - t_old - dt_) < 1.e-4);

  fail = AdvanceStep_i_(0, t_old, t_new, reinit);
  if (fail) return fail;
  CopyStarToPrimary_(tag_current_, tag_next_, tag_current_, tag_next_);
  fail = AdvanceStep_i_(1, t_old, t_new, reinit);
  return fail;
}


void MPCPermafrostSplitFlux::CommitStep(double t_old, double t_new,
        const Tag& tag)
{
  // NOTE: in AJC code, these were flipped.  I believe this is correct, but
  // might be worth checking to see which works better.  Note it should result
  // in larget timestep sizes to get this right, but the physics shouldn't
  // change if this is wrong.  In AJC code, the below comment was still
  // present...
  //
  // Also, if this _did_ need to be flipped, than the call to
  // sub_pks_[i]->CommitStep() in the Advance_Subcycled would be incorrect, and
  // we would have to unravel that loop and not call Commit on the star system.
  //
  // Commit before copy to ensure record for extrapolation in star system uses
  // its own solutions
  MPCSubcycled::CommitStep(t_old, t_new, tag);

  // Copy the primary into the star to advance
  CopyPrimaryToStar_(tag, tag);
}


void
MPCPermafrostSplitFlux::CopyPrimaryToStar_(const Tag& primary,
                                    const Tag& star)
{
  if (is_domain_set_) {
    CopyPrimaryToStar_DomainSet_(primary, star);
  } else {
    CopyPrimaryToStar_Standard_(primary, star);
  }
}


void
MPCPermafrostSplitFlux::CopyStarToPrimary_(
  const Tag& star_current, const Tag& star_next,
  const Tag& primary_current, const Tag& primary_next)
{
  if (is_domain_set_) {
    if (coupling_ == "pressure") CopyStarToPrimary_DomainSet_Pressure_(star_current, star_next, primary_current, primary_next);
    else if (coupling_ == "flux") CopyStarToPrimary_DomainSet_Flux_(star_current, star_next, primary_current, primary_next);
    else if (coupling_ == "hybrid") CopyStarToPrimary_DomainSet_Hybrid_(star_current, star_next, primary_current, primary_next);
    else AMANZI_ASSERT(false);
  } else {
    if (coupling_ == "pressure") CopyStarToPrimary_Standard_Pressure_(star_current, star_next, primary_current, primary_next);
    else if (coupling_ == "flux") CopyStarToPrimary_Standard_Flux_(star_current, star_next, primary_current, primary_next);
    else if (coupling_ == "hybrid") CopyStarToPrimary_Standard_Hybrid_(star_current, star_next, primary_current, primary_next);
    else AMANZI_ASSERT(false);
  }
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system assuming a 3D subsurface domain
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyPrimaryToStar_Standard_(const Tag& primary, const Tag& star)
{
  // copy p primary variable
  auto p_owner = S_->GetRecord(p_primary_variable_star_, star).owner();
  auto& p_star = *S_->GetW<CompositeVector>(p_primary_variable_star_, star, p_owner)
                  .ViewComponent("cell",false);
  const auto& p = *S_->Get<CompositeVector>(p_primary_variable_, primary).ViewComponent("cell",false);
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p[0][c] <= 101325.0) {
      p_star[0][c] = 101325.;
    } else {
      p_star[0][c] = p[0][c];
    }
  }
  ChangedEvaluatorPrimary(p_primary_variable_star_, star, *S_);

  // copy T primary variable
  auto T_owner = S_->GetRecord(T_primary_variable_star_, star).owner();
  auto& T_star = *S_->GetW<CompositeVector>(T_primary_variable_star_, star,
          T_owner).ViewComponent("cell",false);
  const auto& T = *S_->Get<CompositeVector>(T_primary_variable_, primary)
    .ViewComponent("cell",false);
  T_star = T;
  ChangedEvaluatorPrimary(T_primary_variable_star_, star, *S_);
}

// -----------------------------------------------------------------------------
// Copy the primary variable to the star system assuming a DomainSet domain
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyPrimaryToStar_DomainSet_(const Tag& primary, const Tag& star)
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // copy p primary variables into star primary variable
  auto p_owner = S_->GetRecord(p_primary_variable_star_, star).owner();
  auto& p_star = *S_->GetW<CompositeVector>(p_primary_variable_star_, star,
          p_owner).ViewComponent("cell",false);

  auto T_owner = S_->GetRecord(T_primary_variable_star_, star).owner();
  auto& T_star = *S_->GetW<CompositeVector>(T_primary_variable_star_, star,
          T_owner).ViewComponent("cell",false);

  auto ds_iter = domain_set.begin();
  for (int c=0; c!=p_star.MyLength(); ++c) {
    Key p_key = Keys::getKey(*ds_iter, p_primary_variable_suffix_);
    const auto& p = *S_->Get<CompositeVector>(p_key, primary).ViewComponent("cell",false);
    AMANZI_ASSERT(p.MyLength() == 1);
    if (p[0][0] <= 101325.0) {
      p_star[0][c] = 101325.;
    } else {
      p_star[0][c] = p[0][0];
    }

    Key T_key = Keys::getKey(*ds_iter, T_primary_variable_suffix_);
    const auto& T = *S_->Get<CompositeVector>(T_key, primary).ViewComponent("cell",false);
    AMANZI_ASSERT(T.MyLength() == 1);
    T_star[0][c] = T[0][0];
    ++ds_iter;
  }
  ChangedEvaluatorPrimary(p_primary_variable_star_, star, *S_);
  ChangedEvaluatorPrimary(T_primary_variable_star_, star, *S_);
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_Standard_Pressure_(
  const Tag& star_current, const Tag& star_next,
  const Tag& primary_current, const Tag& primary_next)
{
  // copy p primary variables from star primary variable
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, star_next)
    ->ViewComponent("cell", false);
  const auto& WC_star = *S_->GetPtr<CompositeVector>(p_conserved_variable_star_, star_next)
    ->ViewComponent("cell", false);

  auto p_owner = S_->GetRecord(p_primary_variable_, primary_current).owner();
  auto& p = *S_->GetW<CompositeVector>(p_primary_variable_, primary_current,
          p_owner).ViewComponent("cell", false);
  auto& WC = *S_->GetW<CompositeVector>(p_conserved_variable_, primary_current,
          p_conserved_variable_).ViewComponent("cell", false);

  double p_atm_plus_eps = 101325. + 1.e-7;
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > p_atm_plus_eps) {
      p[0][c] = p_star[0][c];
      WC[0][c] = WC_star[0][c];
    }
  }

  // ETC: Probably need to also copy over water content now that we don't have
  // evaluators at the old time?  In this first cut this is a direct
  // translation, so this may not work as is!
  //
  // mark p and subsurface p as changed
  ChangedEvaluatorPrimary(p_primary_variable_, primary_current, *S_);
  // ChangedEvaluatorPrimary(p_conserved_variable_, primary_current, *S_);
  auto p_sub_owner = S_->GetRecord(p_sub_primary_variable_, primary_current).owner();
  CopySurfaceToSubsurface(S_->Get<CompositeVector>(p_primary_variable_, primary_current),
                          S_->GetW<CompositeVector>(p_sub_primary_variable_, primary_current,
                                  p_sub_owner));
  ChangedEvaluatorPrimary(p_sub_primary_variable_, primary_current, *S_);

  // copy T primary variables from star primary variable
  const auto& T_star = *S_->Get<CompositeVector>(T_primary_variable_star_, star_next)
    .ViewComponent("cell", false);
  const auto& E_star = *S_->Get<CompositeVector>(T_conserved_variable_star_, star_next)
    .ViewComponent("cell", false);
  auto T_owner = S_->GetRecord(T_primary_variable_, primary_current).owner();
  auto& T = *S_->GetW<CompositeVector>(T_primary_variable_, primary_current, T_primary_variable_)
    .ViewComponent("cell", false);
  auto& E = *S_->GetW<CompositeVector>(T_conserved_variable_, primary_current, T_conserved_variable_)
    .ViewComponent("cell", false);
  T = T_star;
  E = E_star;

  // mark T and subsurface T as changed
  ChangedEvaluatorPrimary(T_primary_variable_, primary_current, *S_);
  // ChangedEvaluatorPrimary(T_conserved_variable_, primary_current, *S_);
  auto T_sub_owner = S_->GetRecord(T_sub_primary_variable_, primary_current).owner();
  CopySurfaceToSubsurface(S_->Get<CompositeVector>(T_primary_variable_, primary_current),
                          S_->GetW<CompositeVector>(T_sub_primary_variable_, primary_current,
                                  T_sub_owner));
  ChangedEvaluatorPrimary(T_sub_primary_variable_, primary_current, *S_);
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_Standard_Flux_(
  const Tag& star_current, const Tag& star_next,
  const Tag& primary_current, const Tag& primary_next)
{
  double dt = S_->get_time(star_next) - S_->get_time(star_current);

  // mass
  // -- grab the data, difference
  auto& q_div = *S_->GetW<CompositeVector>(p_lateral_flow_source_, primary_next, p_lateral_flow_source_)
                .ViewComponent("cell",false);
  q_div.Update(1.0/dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, star_next).ViewComponent("cell",false),
               -1.0/dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, star_current).ViewComponent("cell",false),
               0.);

  // -- scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_->Get<CompositeVector>(cv_key_, primary_next).ViewComponent("cell",false), q_div, 0.);

  // -- mark the source evaluator as changed to ensure the total source gets updated.
  ChangedEvaluatorPrimary(p_lateral_flow_source_, primary_next, *S_);

  // energy
  // -- grab the data, difference
  auto& qE_div = *S_->GetW<CompositeVector>(T_lateral_flow_source_, primary_next, T_lateral_flow_source_)
                .ViewComponent("cell",false);
  qE_div.Update(1.0/dt,
               *S_->Get<CompositeVector>(T_conserved_variable_star_, star_next).ViewComponent("cell",false),
               -1.0/dt,
               *S_->Get<CompositeVector>(T_conserved_variable_star_, star_current).ViewComponent("cell",false),
               0.);

  // -- scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_->Get<CompositeVector>(cv_key_, primary_next).ViewComponent("cell",false), qE_div, 0.);

  // -- mark the source evaluator as changed to ensure the total source gets updated.
  ChangedEvaluatorPrimary(T_lateral_flow_source_, primary_next, *S_);
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_Standard_Hybrid_(
  const Tag& star_current, const Tag& star_next,
  const Tag& primary_current, const Tag& primary_next)
{
  double dt = S_->get_time(star_next) - S_->get_time(star_current);

  // mass
  // -- grab the data, difference
  auto& q_div = *S_->GetW<CompositeVector>(p_lateral_flow_source_, primary_next, p_lateral_flow_source_)
                .ViewComponent("cell",false);
  q_div.Update(1.0/dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, star_next).ViewComponent("cell",false),
               -1.0/dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, star_current).ViewComponent("cell",false),
               0.);

  // -- scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_->Get<CompositeVector>(cv_key_, primary_next).ViewComponent("cell",false), q_div, 0.);

  // energy
  // -- grab the data, difference
  auto& qE_div = *S_->GetW<CompositeVector>(T_lateral_flow_source_, primary_next, T_lateral_flow_source_)
                .ViewComponent("cell",false);
  qE_div.Update(1.0/dt,
               *S_->Get<CompositeVector>(T_conserved_variable_star_, star_next).ViewComponent("cell",false),
               -1.0/dt,
               *S_->Get<CompositeVector>(T_conserved_variable_star_, star_current).ViewComponent("cell",false),
               0.);

  // -- scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_->Get<CompositeVector>(cv_key_, primary_next).ViewComponent("cell",false), qE_div, 0.);

  // grab the pressure and temp from the star system as well
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, star_next)
    ->ViewComponent("cell", false);
  const auto& T_star = *S_->Get<CompositeVector>(T_primary_variable_star_, star_next)
    .ViewComponent("cell", false);
  const auto& WC_star = *S_->GetPtr<CompositeVector>(p_conserved_variable_star_, star_next)
    ->ViewComponent("cell", false);
  const auto& E_star = *S_->Get<CompositeVector>(T_conserved_variable_star_, star_next)
    .ViewComponent("cell", false);

  // and from the surface system
  auto p_owner = S_->GetRecord(p_primary_variable_, primary_current).owner();
  auto& p = *S_->GetW<CompositeVector>(p_primary_variable_, primary_current,
          p_owner).ViewComponent("cell", false);
  auto& WC = *S_->GetW<CompositeVector>(p_conserved_variable_, primary_current,
          p_conserved_variable_).ViewComponent("cell", false);
  auto T_owner = S_->GetRecord(T_primary_variable_, primary_current).owner();
  auto& T = *S_->GetW<CompositeVector>(T_primary_variable_, primary_current, T_owner)
    .ViewComponent("cell", false);
  auto& E = *S_->GetW<CompositeVector>(T_conserved_variable_, primary_current,
          T_conserved_variable_).ViewComponent("cell", false);

  // in the case of water loss, use pressure.  in the case of water gain, use flux.
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > 101325. && q_div[0][c] < 0.) {
      // use the Dirichlet
      p[0][c] = p_star[0][c];
      WC[0][c] = WC_star[0][c];
      T[0][c] = T_star[0][c];
      E[0][c] = E_star[0][c];

      // set the lateral flux to 0
      q_div[0][c] = 0.;
      qE_div[0][c] = 0.;
    } // otherwise, use the flux, so nothing changes
  }

  // mark both the primary evals and flux evals as changed
  ChangedEvaluatorPrimary(p_lateral_flow_source_, primary_next, *S_);
  ChangedEvaluatorPrimary(T_lateral_flow_source_, primary_next, *S_);

  // mark p and subsurface p as changed
  ChangedEvaluatorPrimary(p_primary_variable_, primary_current, *S_);
  // ChangedEvaluatorPrimary(p_conserved_variable_, primary_current, *S_);
  auto p_sub_owner = S_->GetRecord(p_sub_primary_variable_, primary_current).owner();
  CopySurfaceToSubsurface(S_->Get<CompositeVector>(p_primary_variable_, primary_current),
                          S_->GetW<CompositeVector>(p_sub_primary_variable_, primary_current,
                                  p_sub_owner));
  ChangedEvaluatorPrimary(p_sub_primary_variable_, primary_current, *S_);

  // mark T and subsurface T as changed
  auto T_sub_owner = S_->GetRecord(T_sub_primary_variable_, primary_current).owner();
  ChangedEvaluatorPrimary(T_primary_variable_, primary_current, *S_);
  // ChangedEvaluatorPrimary(T_conserved_variable_, primary_current, *S_);
  CopySurfaceToSubsurface(S_->Get<CompositeVector>(T_primary_variable_, primary_current),
                          S_->GetW<CompositeVector>(T_sub_primary_variable_, primary_current,
                                  T_sub_owner));
  ChangedEvaluatorPrimary(T_sub_primary_variable_, primary_current, *S_);
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system, using the pressure coupling schem
//
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_DomainSet_Pressure_(
  const Tag& star_current, const Tag& star_next,
  const Tag& primary_current, const Tag& primary_next)
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // copy p primary variables into star primary variable
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, star_next)
    ->ViewComponent("cell", false);
  const auto& T_star = *S_->GetPtr<CompositeVector>(T_primary_variable_star_, star_next)
    ->ViewComponent("cell", false);

  auto ds_iter = domain_set.begin();
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > 101325.0000001) {
      Key p_key = Keys::getKey(*ds_iter, p_primary_variable_suffix_);
      auto p_owner = S_->GetRecord(p_key, primary_current).owner();
      auto& p = *S_->GetW<CompositeVector>(p_key, primary_current, p_owner).ViewComponent("cell",false);
      AMANZI_ASSERT(p.MyLength() == 1);
      p[0][0] = p_star[0][c];

      // ?? what about WC?
      ChangedEvaluatorPrimary(p_key, primary_current, *S_);
      Key p_sub_key = Keys::getKey(domain_sub_,
              Keys::getDomainSetIndex(*ds_iter), p_sub_primary_variable_suffix_);
      auto p_sub_owner = S_->GetRecord(p_sub_key, primary_current).owner();
      CopySurfaceToSubsurface(S_->Get<CompositeVector>(p_key, primary_current),
              S_->GetW<CompositeVector>(p_sub_key, primary_current, p_sub_owner));
    }

    Key T_key = Keys::getKey(*ds_iter, T_primary_variable_suffix_);
    auto T_owner = S_->GetRecord(T_key, primary_current).owner();
    auto& T = *S_->GetW<CompositeVector>(T_key, primary_current, T_owner).ViewComponent("cell",false);
    AMANZI_ASSERT(T.MyLength() == 1);
    T[0][0] = T_star[0][c];

    // ?? what about WC?
    ChangedEvaluatorPrimary(T_key, primary_current, *S_);
    Key T_sub_key = Keys::getKey(domain_sub_,
            Keys::getDomainSetIndex(*ds_iter), T_sub_primary_variable_suffix_);
    auto T_sub_owner = S_->GetRecord(T_sub_key, primary_current).owner();
    CopySurfaceToSubsurface(S_->Get<CompositeVector>(T_key, primary_current),
                            S_->GetW<CompositeVector>(T_sub_key, primary_current, T_sub_owner));
    ++ds_iter;
  }
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_DomainSet_Flux_(
  const Tag& star_current, const Tag& star_next,
  const Tag& primary_current, const Tag& primary_next)
{
  double dt = S_->get_time(star_next) - S_->get_time(star_current);
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // grab the data, difference
  Epetra_MultiVector q_div(*S_->Get<CompositeVector>(p_conserved_variable_star_, star_next).ViewComponent("cell",false));
  q_div.Update(1.0/dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, star_next).ViewComponent("cell",false),
               -1.0/dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, star_current).ViewComponent("cell",false),
               0.);
  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_->Get<CompositeVector>(cv_key_, primary_next).ViewComponent("cell",false), q_div, 0.);

  // grab the data, difference
  Epetra_MultiVector qE_div(*S_->Get<CompositeVector>(T_conserved_variable_star_, star_next).ViewComponent("cell",false));
  qE_div.Update(1.0/dt,
               *S_->Get<CompositeVector>(T_conserved_variable_star_, star_next).ViewComponent("cell",false),
               -1.0/dt,
               *S_->Get<CompositeVector>(T_conserved_variable_star_, star_current).ViewComponent("cell",false),
               0.);
  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_->Get<CompositeVector>(cv_key_, primary_next).ViewComponent("cell",false), qE_div, 0.);

  // copy into columns
  auto ds_iter = domain_set.begin();
  for (int c=0; c!=q_div.MyLength(); ++c) {
    Key p_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
    auto p_owner = S_->GetRecord(p_key, primary_next).owner();
    (*S_->GetW<CompositeVector>(p_key, primary_next, p_owner).ViewComponent("cell",false))[0][0] = q_div[0][c];
    ChangedEvaluatorPrimary(p_key, primary_next, *S_);

    Key T_key = Keys::getKey(*ds_iter, T_lateral_flow_source_suffix_);
    auto T_owner = S_->GetRecord(T_key, primary_next).owner();
    (*S_->GetW<CompositeVector>(T_key, primary_next, T_owner).ViewComponent("cell",false))[0][0] = qE_div[0][c];
    ChangedEvaluatorPrimary(T_key, primary_next, *S_);
    ++ds_iter;
  }
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_DomainSet_Hybrid_(
  const Tag& star_current, const Tag& star_next,
  const Tag& primary_current, const Tag& primary_next)
{
  double dt = S_->get_time(star_next) - S_->get_time(star_current);
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // grab the data, difference
  Epetra_MultiVector q_div(*S_->Get<CompositeVector>(p_conserved_variable_star_, star_next).ViewComponent("cell",false));
  q_div.Update(1.0/dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, star_next).ViewComponent("cell",false),
               -1.0/dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, star_current).ViewComponent("cell",false),
               0.);
  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_->Get<CompositeVector>(cv_key_, primary_next).ViewComponent("cell",false), q_div, 0.);

  // grab the data, difference
  Epetra_MultiVector qE_div(*S_->Get<CompositeVector>(T_conserved_variable_star_, star_next).ViewComponent("cell",false));
  qE_div.Update(1.0/dt,
               *S_->Get<CompositeVector>(T_conserved_variable_star_, star_next).ViewComponent("cell",false),
               -1.0/dt,
               *S_->Get<CompositeVector>(T_conserved_variable_star_, star_current).ViewComponent("cell",false),
               0.);
  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_->Get<CompositeVector>(cv_key_, primary_next).ViewComponent("cell",false), qE_div, 0.);

  // copy p primary variables into star primary variable
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, star_next)
    ->ViewComponent("cell", false);
  const auto& T_star = *S_->GetPtr<CompositeVector>(T_primary_variable_star_, star_next)
    ->ViewComponent("cell", false);

  // in the case of water loss, use pressure.  in the case of water gain, use flux.
  auto ds_iter = domain_set.begin();
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > 101325. && q_div[0][c] < 0.) {
      // use the Dirichlet
      Key p_key = Keys::getKey(*ds_iter, p_primary_variable_suffix_);
      auto p_owner = S_->GetRecord(p_key, primary_current).owner();
      auto& p = *S_->GetW<CompositeVector>(p_key, primary_current, p_owner).ViewComponent("cell",false);
      AMANZI_ASSERT(p.MyLength() == 1);
      p[0][0] = p_star[0][c];

      Key T_key = Keys::getKey(*ds_iter, T_primary_variable_suffix_);
      auto T_owner = S_->GetRecord(T_key, primary_current).owner();
      auto& T = *S_->GetW<CompositeVector>(T_key, primary_current, T_owner).ViewComponent("cell",false);
      AMANZI_ASSERT(T.MyLength() == 1);
      T[0][0] = T_star[0][c];

      // ?? what about WC?
      ChangedEvaluatorPrimary(p_key, primary_current, *S_);
      Key p_sub_key = Keys::getKey(domain_sub_,
              Keys::getDomainSetIndex(*ds_iter), p_sub_primary_variable_suffix_);
      auto p_sub_owner = S_->GetRecord(p_sub_key, primary_current).owner();
      CopySurfaceToSubsurface(S_->Get<CompositeVector>(p_key, primary_current),
              S_->GetW<CompositeVector>(p_sub_key, primary_current, p_sub_owner));

      // ?? what about WC?
      ChangedEvaluatorPrimary(T_key, primary_current, *S_);
      Key T_sub_key = Keys::getKey(domain_sub_,
              Keys::getDomainSetIndex(*ds_iter), T_sub_primary_variable_suffix_);
      auto T_sub_owner = S_->GetRecord(T_sub_key, primary_current).owner();
      CopySurfaceToSubsurface(S_->Get<CompositeVector>(T_key, primary_current),
              S_->GetW<CompositeVector>(T_sub_key, primary_current, T_sub_owner));

      // set the lateral flux to 0
      Key p_lf_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
      (*S_->GetW<CompositeVector>(p_lf_key, primary_next, p_lf_key).ViewComponent("cell",false))[0][0] = 0.;
      ChangedEvaluatorPrimary(p_lf_key, primary_next, *S_);

      Key T_lf_key = Keys::getKey(*ds_iter, T_lateral_flow_source_suffix_);
      (*S_->GetW<CompositeVector>(T_lf_key, primary_next, T_lf_key).ViewComponent("cell",false))[0][0] = 0.;
      ChangedEvaluatorPrimary(T_lf_key, primary_next, *S_);

    } else {
      // use flux
      Key p_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
      auto p_owner = S_->GetRecord(p_key, primary_next).owner();
      (*S_->GetW<CompositeVector>(p_key, primary_next, p_owner).ViewComponent("cell",false))[0][0] = q_div[0][c];
      ChangedEvaluatorPrimary(p_key, primary_next, *S_);

      Key T_key = Keys::getKey(*ds_iter, T_lateral_flow_source_suffix_);
      auto T_owner = S_->GetRecord(T_key, primary_next).owner();
      (*S_->GetW<CompositeVector>(T_key, primary_next, T_owner).ViewComponent("cell",false))[0][0] = qE_div[0][c];
      ChangedEvaluatorPrimary(T_key, primary_next, *S_);
    }
    ++ds_iter;
  }
}


} // namespace


