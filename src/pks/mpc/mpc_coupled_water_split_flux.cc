/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "mpc_coupled_water_split_flux.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "pk_helpers.hh"

namespace Amanzi {

MPCCoupledWaterSplitFlux::MPCCoupledWaterSplitFlux(
  Teuchos::ParameterList& FElist,
  const Teuchos::RCP<Teuchos::ParameterList>& plist,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& solution)
  : PK(FElist, plist, S, solution), MPCSubcycled(FElist, plist, S, solution)
{
  // collect domain names
  domain_set_ = Keys::readDomain(*plist_);          // e.g. surface or surface_column:*
  domain_star_ = Keys::readDomain(*plist_, "star"); // e.g. surface_star

  // determine whether we are coupling subdomains or coupling 3D domains
  is_domain_set_ = S_->HasDomainSet(domain_set_);
  if (is_domain_set_)
    domain_ = Keys::getDomainInSet(domain_set_, "*");
  else
    domain_ = domain_set_;

  domain_sub_ = Keys::readDomainHint(*plist_, domain_set_, "surface", "subsurface");

  // determine the coupling strategy: "pressure" passes the pressure field,
  // "flux" the flux field, while "hybrid" passes one or the other depending
  // upon conditions.  "hybrid" is the most robust.
  coupling_ = plist_->get<std::string>("coupling type", "hybrid");
  if (coupling_ != "pressure" && coupling_ != "flux" && coupling_ != "hybrid") {
    Errors::Message msg("WeakMPCSemiCoupled: \"coupling type\" must be one of \"pressure\", "
                        "\"flux\", or \"hybrid\".");
    Exceptions::amanzi_throw(msg);
  }

  // collect keys and names
  // -- primary variable for the main surface domain
  p_primary_variable_ = Keys::readKey(*plist_, domain_, "pressure primary variable", "pressure");
  p_primary_variable_suffix_ = Keys::getVarName(p_primary_variable_);

  // -- primary variable for the main subsurface domain
  p_sub_primary_variable_ =
    Keys::readKey(*plist_, domain_sub_, "subsurface pressure primary variable", "pressure");
  p_sub_primary_variable_suffix_ = Keys::getVarName(p_sub_primary_variable_);

  // -- need to save the updated conserved quantity too
  p_conserved_variable_ =
    Keys::readKey(*plist_, domain_, "water conserved quantity", "water_content");
  p_conserved_variable_star_ =
    Keys::readKey(*plist_, domain_star_, "water conserved quantity star", "water_content");

  // -- primary variable for the star domain
  p_primary_variable_star_ = Keys::readKey(
    *plist_, domain_star_, "pressure primary variable star", Keys::getVarName(p_primary_variable_));

  // -- flux variables for coupling
  if (coupling_ != "pressure") {
    p_lateral_flow_source_ =
      Keys::readKey(*plist_, domain_, "water lateral flow source", "water_lateral_flow_source");
    p_lateral_flow_source_suffix_ = Keys::getVarName(p_lateral_flow_source_);
    cv_key_ = Keys::readKey(*plist_, domain_star_, "cell volume", "cell_volume");
  }
};


void
MPCCoupledWaterSplitFlux::Setup()
{
  MPCSubcycled::Setup();

  if (coupling_ != "pressure") {
    if (is_domain_set_) {
      auto domain_set = S_->GetDomainSet(domain_set_);
      for (const auto& domain : *domain_set) {
        auto p_key = Keys::getKey(domain, p_lateral_flow_source_suffix_);
        Tag ds_tag_next = get_ds_tag_next_(domain);
        S_->Require<CompositeVector, CompositeVectorSpace>(p_key, ds_tag_next, p_key)
          .SetMesh(S_->GetMesh(domain))
          ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
        requireEvaluatorPrimary(p_key, ds_tag_next, *S_);
      }
    } else {
      S_->Require<CompositeVector, CompositeVectorSpace>(
          p_lateral_flow_source_, tags_[1].second, p_lateral_flow_source_)
        .SetMesh(S_->GetMesh(Keys::getDomain(p_lateral_flow_source_)))
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      requireEvaluatorPrimary(p_lateral_flow_source_, tags_[1].second, *S_);
    }

    // also need conserved quantities at old and new times
    S_->Require<CompositeVector, CompositeVectorSpace>(p_conserved_variable_star_, tags_[0].second)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->Require<CompositeVector, CompositeVectorSpace>(p_conserved_variable_star_, tags_[0].first)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    S_->RequireEvaluator(p_conserved_variable_star_, tags_[0].second);
    //S_->RequireEvaluator(p_conserved_variable_star_, tags_[0].first);
  }
}


void
MPCCoupledWaterSplitFlux::Initialize()
{
  sub_pks_[1]->Initialize();
  CopyPrimaryToStar_();
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
  auto p_owner = S_->GetRecord(p_primary_variable_star_, tags_[0].second).owner();
  S_->GetRecordW(p_primary_variable_star_, tags_[0].second, p_owner).set_initialized();

  auto& pstar = *S_->Get<CompositeVector>(p_primary_variable_star_, tags_[0].second)
                   .ViewComponent("cell", false);

  // set the fluxes as initialized -- they will get set by later calls to
  // CopyStarToPrimary, so no need for values, but do need to toggle the flag.
  if (coupling_ != "pressure") {
    if (is_domain_set_) {
      auto domain_set = S_->GetDomainSet(domain_set_);
      for (const auto& domain : *domain_set) {
        auto pkey = Keys::getKey(domain, p_lateral_flow_source_suffix_);
        Tag ds_tag_next = get_ds_tag_next_(domain);
        auto p_owner = S_->GetRecord(pkey, ds_tag_next).owner();
        S_->GetRecordW(pkey, ds_tag_next, p_owner).set_initialized();
      }
    } else {
      S_->GetRecordW(p_lateral_flow_source_, tags_[1].second, p_lateral_flow_source_)
        .set_initialized();
    }
  }

  int i = 0;
  for (const auto& tag : tags_) {
    if (subcycling_[i]) S_->GetRecordW("dt", tag.second, name()).set_initialized();
    ++i;
  }
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool
MPCCoupledWaterSplitFlux::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  fail = AdvanceStep_i_(0, t_old, t_new, reinit);
  if (fail) return fail;
  CopyStarToPrimary_();
  fail = AdvanceStep_i_(1, t_old, t_new, reinit);
  return fail;
}


void
MPCCoupledWaterSplitFlux::CommitStep(double t_old, double t_new, const Tag& tag)
{
  // Copy the primary into the star to advance.
  //
  // Note that this must be done prior to commit so that all copies of the star
  // system get the right values to start the new timestep.
  CopyPrimaryToStar_();

  MPCSubcycled::CommitStep(t_old, t_new, tag);
}


void
MPCCoupledWaterSplitFlux::CopyPrimaryToStar_()
{
  if (is_domain_set_) {
    CopyPrimaryToStar_DomainSet_();
  } else {
    CopyPrimaryToStar_Standard_();
  }
}


void
MPCCoupledWaterSplitFlux::CopyStarToPrimary_()
{
  if (is_domain_set_) {
    if (coupling_ == "pressure")
      CopyStarToPrimary_DomainSet_Pressure_();
    else if (coupling_ == "flux")
      CopyStarToPrimary_DomainSet_Flux_();
    else if (coupling_ == "hybrid")
      CopyStarToPrimary_DomainSet_Hybrid_();
    else
      AMANZI_ASSERT(false);
  } else {
    if (coupling_ == "pressure")
      CopyStarToPrimary_Standard_Pressure_();
    else if (coupling_ == "flux")
      CopyStarToPrimary_Standard_Flux_();
    else if (coupling_ == "hybrid")
      CopyStarToPrimary_Standard_Hybrid_();
    else
      AMANZI_ASSERT(false);
  }
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system assuming a 3D subsurface domain
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyPrimaryToStar_Standard_()
{
  // copy p primary variable
  auto p_owner = S_->GetRecord(p_primary_variable_star_, tags_[0].second).owner();
  auto& p_star = *S_->GetW<CompositeVector>(p_primary_variable_star_, tags_[0].second, p_owner)
                    .ViewComponent("cell", false);
  const auto& p =
    *S_->Get<CompositeVector>(p_primary_variable_, tags_[1].second).ViewComponent("cell", false);
  for (int c = 0; c != p_star.MyLength(); ++c) {
    if (p[0][c] <= 101325.0) {
      p_star[0][c] = 101325.;
    } else {
      p_star[0][c] = p[0][c];
    }
  }
  changedEvaluatorPrimary(p_primary_variable_star_, tags_[0].second, *S_);
}

// -----------------------------------------------------------------------------
// Copy the primary variable to the star system assuming a DomainSet domain
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyPrimaryToStar_DomainSet_()
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // copy p primary variables into star primary variable
  auto p_owner = S_->GetRecord(p_primary_variable_star_, tags_[0].second).owner();
  auto& p_star = *S_->GetW<CompositeVector>(p_primary_variable_star_, tags_[0].second, p_owner)
                    .ViewComponent("cell", false);

  auto ds_iter = domain_set.begin();
  for (int c = 0; c != p_star.MyLength(); ++c) {
    Key p_key = Keys::getKey(*ds_iter, p_primary_variable_suffix_);
    Tag ds_tag_next = get_ds_tag_next_(*ds_iter);
    const auto& p = *S_->Get<CompositeVector>(p_key, ds_tag_next).ViewComponent("cell", false);
    AMANZI_ASSERT(p.MyLength() == 1);
    if (p[0][0] <= 101325.0) {
      p_star[0][c] = 101325.;
    } else {
      p_star[0][c] = p[0][0];
    }
    ++ds_iter;
  }
  changedEvaluatorPrimary(p_primary_variable_star_, tags_[0].second, *S_);
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyStarToPrimary_Standard_Pressure_()
{
  // copy p primary variables from star primary variable
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, tags_[0].second)
                          ->ViewComponent("cell", false);
  const auto& WC_star = *S_->GetPtr<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                           ->ViewComponent("cell", false);

  auto p_owner = S_->GetRecord(p_primary_variable_, tags_[1].first).owner();
  auto& p = *S_->GetW<CompositeVector>(p_primary_variable_, tags_[1].first, p_owner)
               .ViewComponent("cell", false);
  auto& WC =
    *S_->GetW<CompositeVector>(p_conserved_variable_, tags_[1].first, p_conserved_variable_)
       .ViewComponent("cell", false);

  double p_atm_plus_eps = 101325. + 1.e-7;
  for (int c = 0; c != p_star.MyLength(); ++c) {
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
  changedEvaluatorPrimary(p_primary_variable_, tags_[1].first, *S_);
  // changedEvaluatorPrimary(p_conserved_variable_, tags_[1].first, *S_);
  auto p_sub_owner = S_->GetRecord(p_sub_primary_variable_, tags_[1].first).owner();
  CopySurfaceToSubsurface(
    S_->Get<CompositeVector>(p_primary_variable_, tags_[1].first),
    S_->GetW<CompositeVector>(p_sub_primary_variable_, tags_[1].first, p_sub_owner));
  changedEvaluatorPrimary(p_sub_primary_variable_, tags_[1].first, *S_);
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyStarToPrimary_Standard_Flux_()
{
  double dt = S_->get_time(tags_[0].second) - S_->get_time(tags_[0].first);

  // mass
  // -- grab the data, difference
  auto& q_div =
    *S_->GetW<CompositeVector>(p_lateral_flow_source_, tags_[1].second, p_lateral_flow_source_)
       .ViewComponent("cell", false);
  q_div.Update(1.0 / dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                  .ViewComponent("cell", false),
               -1.0 / dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].first)
                  .ViewComponent("cell", false),
               0.);

  // -- scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[1].second).ViewComponent("cell", false),
    q_div,
    0.);

  // -- mark the source evaluator as changed to ensure the total source gets updated.
  changedEvaluatorPrimary(p_lateral_flow_source_, tags_[1].second, *S_);
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyStarToPrimary_Standard_Hybrid_()
{
  double dt = S_->get_time(tags_[0].second) - S_->get_time(tags_[0].first);

  // mass
  // -- grab the data, difference
  auto& q_div =
    *S_->GetW<CompositeVector>(p_lateral_flow_source_, tags_[1].second, p_lateral_flow_source_)
       .ViewComponent("cell", false);
  q_div.Update(1.0 / dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                  .ViewComponent("cell", false),
               -1.0 / dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].first)
                  .ViewComponent("cell", false),
               0.);

  // -- scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[1].second).ViewComponent("cell", false),
    q_div,
    0.);

  // grab the pressure and temp from the star system as well
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, tags_[0].second)
                          ->ViewComponent("cell", false);
  const auto& WC_star = *S_->GetPtr<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                           ->ViewComponent("cell", false);

  // and from the surface system
  auto p_owner = S_->GetRecord(p_primary_variable_, tags_[1].first).owner();
  auto& p = *S_->GetW<CompositeVector>(p_primary_variable_, tags_[1].first, p_owner)
               .ViewComponent("cell", false);
  auto& WC =
    *S_->GetW<CompositeVector>(p_conserved_variable_, tags_[1].first, p_conserved_variable_)
       .ViewComponent("cell", false);

  // in the case of water loss, use pressure.  in the case of water gain, use flux.
  for (int c = 0; c != p_star.MyLength(); ++c) {
    if (p_star[0][c] > 101325. && q_div[0][c] < 0.) {
      // use the Dirichlet
      p[0][c] = p_star[0][c];
      WC[0][c] = WC_star[0][c];

      // set the lateral flux to 0
      q_div[0][c] = 0.;
    } // otherwise, use the flux, so nothing changes
  }

  // mark both the primary evals and flux evals as changed
  changedEvaluatorPrimary(p_lateral_flow_source_, tags_[1].second, *S_);

  // mark p and subsurface p as changed
  changedEvaluatorPrimary(p_primary_variable_, tags_[1].first, *S_);
  // changedEvaluatorPrimary(p_conserved_variable_, tags_[1].first, *S_);
  auto p_sub_owner = S_->GetRecord(p_sub_primary_variable_, tags_[1].first).owner();
  CopySurfaceToSubsurface(
    S_->Get<CompositeVector>(p_primary_variable_, tags_[1].first),
    S_->GetW<CompositeVector>(p_sub_primary_variable_, tags_[1].first, p_sub_owner));
  changedEvaluatorPrimary(p_sub_primary_variable_, tags_[1].first, *S_);
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system, using the pressure coupling schem
//
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyStarToPrimary_DomainSet_Pressure_()
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // copy p primary variables into star primary variable
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, tags_[0].second)
                          ->ViewComponent("cell", false);

  auto ds_iter = domain_set.begin();
  for (int c = 0; c != p_star.MyLength(); ++c) {
    Tag ds_tag_next = get_ds_tag_next_(*ds_iter);
    Tag ds_tag_current = get_ds_tag_current_(*ds_iter);

    if (p_star[0][c] > 101325.0000001) {
      Key p_key = Keys::getKey(*ds_iter, p_primary_variable_suffix_);
      auto p_owner = S_->GetRecord(p_key, ds_tag_current).owner();
      auto& p =
        *S_->GetW<CompositeVector>(p_key, ds_tag_current, p_owner).ViewComponent("cell", false);
      AMANZI_ASSERT(p.MyLength() == 1);
      p[0][0] = p_star[0][c];

      // ?? what about WC?
      changedEvaluatorPrimary(p_key, ds_tag_current, *S_);
      Key p_sub_key = Keys::getKey(
        domain_sub_, Keys::getDomainSetIndex(*ds_iter), p_sub_primary_variable_suffix_);
      auto p_sub_owner = S_->GetRecord(p_sub_key, ds_tag_current).owner();
      CopySurfaceToSubsurface(S_->Get<CompositeVector>(p_key, ds_tag_current),
                              S_->GetW<CompositeVector>(p_sub_key, ds_tag_current, p_sub_owner));
    }
    ++ds_iter;
  }
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyStarToPrimary_DomainSet_Flux_()
{
  double dt = S_->get_time(tags_[0].second) - S_->get_time(tags_[0].first);
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // grab the data, difference
  Epetra_MultiVector q_div(*S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                              .ViewComponent("cell", false));
  q_div.Update(1.0 / dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                  .ViewComponent("cell", false),
               -1.0 / dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].first)
                  .ViewComponent("cell", false),
               0.);
  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[0].second).ViewComponent("cell", false),
    q_div,
    0.);

  // copy into columns
  auto ds_iter = domain_set.begin();
  for (int c = 0; c != q_div.MyLength(); ++c) {
    Tag ds_tag_next = get_ds_tag_next_(*ds_iter);
    Key p_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
    auto p_owner = S_->GetRecord(p_key, ds_tag_next).owner();
    (*S_->GetW<CompositeVector>(p_key, ds_tag_next, p_owner).ViewComponent("cell", false))[0][0] =
      q_div[0][c];
    changedEvaluatorPrimary(p_key, ds_tag_next, *S_);
    ++ds_iter;
  }
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCCoupledWaterSplitFlux::CopyStarToPrimary_DomainSet_Hybrid_()
{
  double dt = S_->get_time(tags_[0].second) - S_->get_time(tags_[0].first);
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // grab the data, difference
  Epetra_MultiVector q_div(*S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                              .ViewComponent("cell", false));
  q_div.Update(1.0 / dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                  .ViewComponent("cell", false),
               -1.0 / dt,
               *S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].first)
                  .ViewComponent("cell", false),
               0.);
  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[0].second).ViewComponent("cell", false),
    q_div,
    0.);

  // copy p primary variables into star primary variable
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, tags_[0].second)
                          ->ViewComponent("cell", false);

  // in the case of water loss, use pressure.  in the case of water gain, use flux.
  auto ds_iter = domain_set.begin();
  for (int c = 0; c != p_star.MyLength(); ++c) {
    Tag ds_tag_next = get_ds_tag_next_(*ds_iter);
    Tag ds_tag_current = get_ds_tag_current_(*ds_iter);
    if (p_star[0][c] > 101325. && q_div[0][c] < 0.) {
      // use the Dirichlet
      Key p_key = Keys::getKey(*ds_iter, p_primary_variable_suffix_);
      auto p_owner = S_->GetRecord(p_key, ds_tag_current).owner();
      auto& p =
        *S_->GetW<CompositeVector>(p_key, ds_tag_current, p_owner).ViewComponent("cell", false);
      AMANZI_ASSERT(p.MyLength() == 1);
      p[0][0] = p_star[0][c];

      // ?? what about WC?
      changedEvaluatorPrimary(p_key, ds_tag_current, *S_);
      Key p_sub_key = Keys::getKey(
        domain_sub_, Keys::getDomainSetIndex(*ds_iter), p_sub_primary_variable_suffix_);
      auto p_sub_owner = S_->GetRecord(p_sub_key, ds_tag_current).owner();
      CopySurfaceToSubsurface(S_->Get<CompositeVector>(p_key, ds_tag_current),
                              S_->GetW<CompositeVector>(p_sub_key, ds_tag_current, p_sub_owner));

      // set the lateral flux to 0
      Key p_lf_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
      (*S_->GetW<CompositeVector>(p_lf_key, ds_tag_next, p_lf_key)
          .ViewComponent("cell", false))[0][0] = 0.;
      changedEvaluatorPrimary(p_lf_key, ds_tag_next, *S_);

    } else {
      // use flux
      Key p_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
      auto p_owner = S_->GetRecord(p_key, ds_tag_next).owner();
      (*S_->GetW<CompositeVector>(p_key, ds_tag_next, p_owner).ViewComponent("cell", false))[0][0] =
        q_div[0][c];
      changedEvaluatorPrimary(p_key, ds_tag_next, *S_);
    }
    ++ds_iter;
  }
}


} // namespace Amanzi
