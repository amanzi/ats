/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "mpc_permafrost_split_flux.hh"
#include "mpc_surface_subsurface_helpers.hh"
#include "pk_helpers.hh"

namespace Amanzi {

MPCPermafrostSplitFlux::MPCPermafrostSplitFlux(Teuchos::ParameterList& FElist,
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
  domain_snow_ = Keys::readDomainHint(*plist_, domain_set_, "surface", "snow");

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
  T_primary_variable_ =
    Keys::readKey(*plist_, domain_, "temperature primary variable", "temperature");
  p_primary_variable_suffix_ = Keys::getVarName(p_primary_variable_);
  T_primary_variable_suffix_ = Keys::getVarName(T_primary_variable_);

  // -- primary variable for the main subsurface domain
  p_sub_primary_variable_ =
    Keys::readKey(*plist_, domain_sub_, "subsurface pressure primary variable", "pressure");
  T_sub_primary_variable_ =
    Keys::readKey(*plist_, domain_sub_, "subsurface temperature primary variable", "temperature");
  p_sub_primary_variable_suffix_ = Keys::getVarName(p_sub_primary_variable_);
  T_sub_primary_variable_suffix_ = Keys::getVarName(T_sub_primary_variable_);

  // -- need to save the updated conserved quantity too
  p_conserved_variable_ =
    Keys::readKey(*plist_, domain_, "water conserved quantity", "water_content");
  p_conserved_variable_suffix_ = Keys::getVarName(p_conserved_variable_);
  p_conserved_variable_star_ =
    Keys::readKey(*plist_, domain_star_, "water conserved quantity star", "water_content");

  T_conserved_variable_ = Keys::readKey(*plist_, domain_, "energy conserved quantity", "energy");
  T_conserved_variable_suffix_ = Keys::getVarName(T_conserved_variable_);
  T_conserved_variable_star_ =
    Keys::readKey(*plist_, domain_star_, "energy conserved quantity star", "energy");

  // -- primary variable for the star domain
  p_primary_variable_star_ = Keys::readKey(
    *plist_, domain_star_, "pressure primary variable star", Keys::getVarName(p_primary_variable_));
  T_primary_variable_star_ = Keys::readKey(*plist_,
                                           domain_star_,
                                           "temperature primary variable star",
                                           Keys::getVarName(T_primary_variable_));

  // -- flux variables for coupling
  if (coupling_ != "pressure") {
    p_lateral_flow_source_ =
      Keys::readKey(*plist_, domain_, "water lateral flow source", "water_lateral_flow_source");
    p_lateral_flow_source_suffix_ = Keys::getVarName(p_lateral_flow_source_);
    T_lateral_flow_source_ =
      Keys::readKey(*plist_, domain_, "energy lateral flow source", "energy_lateral_flow_source");
    T_lateral_flow_source_suffix_ = Keys::getVarName(T_lateral_flow_source_);

    cv_key_ = Keys::readKey(*plist_, domain_star_, "cell volume", "cell_volume");
  }
};


void
MPCPermafrostSplitFlux::Setup()
{
  MPCSubcycled::Setup();

  if (coupling_ != "pressure") {
    if (is_domain_set_) {
      auto domain_set = S_->GetDomainSet(domain_set_);
      for (const auto& domain : *domain_set) {
        auto p_key = Keys::getKey(domain, p_lateral_flow_source_suffix_);
        Tag ds_tag_next = get_ds_tag_next_(domain);
        requireAtNext(p_key, ds_tag_next, *S_, name_)
          .SetMesh(S_->GetMesh(domain))
          ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

        auto T_key = Keys::getKey(domain, T_lateral_flow_source_suffix_);
        requireAtNext(T_key, ds_tag_next, *S_, name_)
          .SetMesh(S_->GetMesh(domain))
          ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
      }
    } else {
      requireAtNext(p_lateral_flow_source_, tags_[1].second, *S_, name_)
        .SetMesh(S_->GetMesh(Keys::getDomain(p_lateral_flow_source_)))
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

      requireAtNext(T_lateral_flow_source_, tags_[1].second, *S_, name_)
        .SetMesh(S_->GetMesh(Keys::getDomain(T_lateral_flow_source_)))
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }

    // also need conserved quantities at old and new times
    requireAtNext(p_conserved_variable_star_, tags_[0].second, *S_)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    requireAtCurrent(p_conserved_variable_star_, tags_[0].first, *S_)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

    requireAtNext(T_conserved_variable_star_, tags_[0].second, *S_)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    requireAtCurrent(T_conserved_variable_star_, tags_[0].first, *S_)
      .SetMesh(S_->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  }
}


void
MPCPermafrostSplitFlux::Initialize()
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
  auto T_owner = S_->GetRecord(T_primary_variable_star_, tags_[0].second).owner();
  S_->GetRecordW(T_primary_variable_star_, tags_[0].second, T_owner).set_initialized();

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
        S_->GetRecordW(pkey, ds_tag_next, name_).set_initialized();

        auto Tkey = Keys::getKey(domain, T_lateral_flow_source_suffix_);
        S_->GetRecordW(Tkey, ds_tag_next, name_).set_initialized();
      }
    } else {
      S_->GetRecordW(p_lateral_flow_source_, tags_[1].second, name_).set_initialized();
      S_->GetRecordW(T_lateral_flow_source_, tags_[1].second, name_).set_initialized();
    }
  }

  int i = 0;
  for (const auto& tag : tags_) {
    if (subcycling_[i]) S_->GetRecordW("dt", tag.first, name()).set_initialized();
    ++i;
  }
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool
MPCPermafrostSplitFlux::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  fail = AdvanceStep_i_(0, t_old, t_new, reinit);
  if (fail) return fail;
  CopyStarToPrimary_();
  fail = AdvanceStep_i_(1, t_old, t_new, reinit);
  return fail;
}


void
MPCPermafrostSplitFlux::CommitStep(double t_old, double t_new, const Tag& tag)
{
  // Copy the primary into the star to advance
  CopyPrimaryToStar_();

  // Note we must copy before we advance to ensure that the star system can
  // move the new pressure into the old pressure's value for that substep.
  //
  // This commits to the star system's tags.
  if (subcycling_[0]) {
    sub_pks_[0]->CommitStep(
      S_->get_time(tags_[0].first), S_->get_time(tags_[0].second), tags_[0].second);
  }

  // finally we can do the global commit
  MPCSubcycled::CommitStep(t_old, t_new, tag);
}


void
MPCPermafrostSplitFlux::CopyPrimaryToStar_()
{
  if (is_domain_set_) {
    CopyPrimaryToStar_DomainSet_();
  } else {
    CopyPrimaryToStar_Standard_();
  }
}


void
MPCPermafrostSplitFlux::CopyStarToPrimary_()
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
MPCPermafrostSplitFlux::CopyPrimaryToStar_Standard_()
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

  // copy T primary variable
  auto T_owner = S_->GetRecord(T_primary_variable_star_, tags_[0].second).owner();
  auto& T_star = *S_->GetW<CompositeVector>(T_primary_variable_star_, tags_[0].second, T_owner)
                    .ViewComponent("cell", false);
  const auto& T =
    *S_->Get<CompositeVector>(T_primary_variable_, tags_[1].second).ViewComponent("cell", false);
  T_star = T;
  changedEvaluatorPrimary(T_primary_variable_star_, tags_[0].second, *S_);
}

// -----------------------------------------------------------------------------
// Copy the primary variable to the star system assuming a DomainSet domain
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyPrimaryToStar_DomainSet_()
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // copy p primary variables into star primary variable
  auto p_owner = S_->GetRecord(p_primary_variable_star_, tags_[0].second).owner();
  auto& p_star = *S_->GetW<CompositeVector>(p_primary_variable_star_, tags_[0].second, p_owner)
                    .ViewComponent("cell", false);

  auto T_owner = S_->GetRecord(T_primary_variable_star_, tags_[0].second).owner();
  auto& T_star = *S_->GetW<CompositeVector>(T_primary_variable_star_, tags_[0].second, T_owner)
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

    Key T_key = Keys::getKey(*ds_iter, T_primary_variable_suffix_);
    const auto& T = *S_->Get<CompositeVector>(T_key, ds_tag_next).ViewComponent("cell", false);
    AMANZI_ASSERT(T.MyLength() == 1);
    T_star[0][c] = T[0][0];
    ++ds_iter;
  }
  changedEvaluatorPrimary(p_primary_variable_star_, tags_[0].second, *S_);
  changedEvaluatorPrimary(T_primary_variable_star_, tags_[0].second, *S_);
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_Standard_Pressure_()
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

  // mark p and WC as changed
  changedEvaluatorPrimary(p_primary_variable_, tags_[1].first, *S_);
  changedEvaluatorPrimary(p_conserved_variable_, tags_[1].first, *S_);

  // copy to subsurface and mark that as changed
  auto p_sub_owner = S_->GetRecord(p_sub_primary_variable_, tags_[1].first).owner();
  CopySurfaceToSubsurface(
    S_->Get<CompositeVector>(p_primary_variable_, tags_[1].first),
    S_->GetW<CompositeVector>(p_sub_primary_variable_, tags_[1].first, p_sub_owner));
  changedEvaluatorPrimary(p_sub_primary_variable_, tags_[1].first, *S_);

  // copy T primary variables from star primary variable
  const auto& T_star = *S_->Get<CompositeVector>(T_primary_variable_star_, tags_[0].second)
                          .ViewComponent("cell", false);
  const auto& E_star = *S_->Get<CompositeVector>(T_conserved_variable_star_, tags_[0].second)
                          .ViewComponent("cell", false);
  auto T_owner = S_->GetRecord(T_primary_variable_, tags_[1].first).owner();
  auto& T = *S_->GetW<CompositeVector>(T_primary_variable_, tags_[1].first, T_primary_variable_)
               .ViewComponent("cell", false);
  auto& E = *S_->GetW<CompositeVector>(T_conserved_variable_, tags_[1].first, T_conserved_variable_)
               .ViewComponent("cell", false);
  T = T_star;
  E = E_star;

  // mark T and E as changed
  changedEvaluatorPrimary(T_primary_variable_, tags_[1].first, *S_);
  changedEvaluatorPrimary(T_conserved_variable_, tags_[1].first, *S_);

  // copy to subsurface and mark that as changed
  auto T_sub_owner = S_->GetRecord(T_sub_primary_variable_, tags_[1].first).owner();
  CopySurfaceToSubsurface(
    S_->Get<CompositeVector>(T_primary_variable_, tags_[1].first),
    S_->GetW<CompositeVector>(T_sub_primary_variable_, tags_[1].first, T_sub_owner));
  changedEvaluatorPrimary(T_sub_primary_variable_, tags_[1].first, *S_);
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_Standard_Flux_()
{
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  // mass
  // -- grab the data, difference
  auto& q_div = *S_->GetW<CompositeVector>(p_lateral_flow_source_, tags_[1].second, name_)
                   .ViewComponent("cell", false);
  q_div.Update(
    1.0 / dt,
    *S_->Get<CompositeVector>(p_conserved_variable_star_, tag_next_).ViewComponent("cell", false),
    -1.0 / dt,
    *S_->Get<CompositeVector>(p_conserved_variable_star_, tag_current_)
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

  // energy
  // -- grab the data, difference
  auto& qE_div = *S_->GetW<CompositeVector>(T_lateral_flow_source_, tags_[1].second, name_)
                    .ViewComponent("cell", false);
  qE_div.Update(
    1.0 / dt,
    *S_->Get<CompositeVector>(T_conserved_variable_star_, tag_next_).ViewComponent("cell", false),
    -1.0 / dt,
    *S_->Get<CompositeVector>(T_conserved_variable_star_, tag_current_)
       .ViewComponent("cell", false),
    0.);

  // -- scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[1].second).ViewComponent("cell", false),
    qE_div,
    0.);

  // -- mark the source evaluator as changed to ensure the total source gets updated.
  changedEvaluatorPrimary(T_lateral_flow_source_, tags_[1].second, *S_);
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_Standard_Hybrid_()
{
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  // mass
  // -- grab the data, difference
  auto& q_div = *S_->GetW<CompositeVector>(p_lateral_flow_source_, tags_[1].second, name_)
                   .ViewComponent("cell", false);
  q_div.Update(
    1.0 / dt,
    *S_->Get<CompositeVector>(p_conserved_variable_star_, tag_next_).ViewComponent("cell", false),
    -1.0 / dt,
    *S_->Get<CompositeVector>(p_conserved_variable_star_, tag_current_)
       .ViewComponent("cell", false),
    0.);

  // -- scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[1].second).ViewComponent("cell", false),
    q_div,
    0.);

  // energy
  // -- grab the data, difference
  auto& qE_div = *S_->GetW<CompositeVector>(T_lateral_flow_source_, tags_[1].second, name_)
                    .ViewComponent("cell", false);
  qE_div.Update(
    1.0 / dt,
    *S_->Get<CompositeVector>(T_conserved_variable_star_, tag_next_).ViewComponent("cell", false),
    -1.0 / dt,
    *S_->Get<CompositeVector>(T_conserved_variable_star_, tag_current_)
       .ViewComponent("cell", false),
    0.);

  // -- scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[1].second).ViewComponent("cell", false),
    qE_div,
    0.);

  // grab the pressure and temp from the star system as well
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, tags_[0].second)
                          ->ViewComponent("cell", false);
  const auto& T_star = *S_->Get<CompositeVector>(T_primary_variable_star_, tags_[0].second)
                          .ViewComponent("cell", false);
  const auto& WC_star = *S_->GetPtr<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                           ->ViewComponent("cell", false);
  const auto& E_star = *S_->Get<CompositeVector>(T_conserved_variable_star_, tags_[0].second)
                          .ViewComponent("cell", false);

  // and from the surface system
  auto p_owner = S_->GetRecord(p_primary_variable_, tags_[1].first).owner();
  auto& p = *S_->GetW<CompositeVector>(p_primary_variable_, tags_[1].first, p_owner)
               .ViewComponent("cell", false);
  auto& WC =
    *S_->GetW<CompositeVector>(p_conserved_variable_, tags_[1].first, p_conserved_variable_)
       .ViewComponent("cell", false);
  auto T_owner = S_->GetRecord(T_primary_variable_, tags_[1].first).owner();
  auto& T = *S_->GetW<CompositeVector>(T_primary_variable_, tags_[1].first, T_owner)
               .ViewComponent("cell", false);
  auto& E = *S_->GetW<CompositeVector>(T_conserved_variable_, tags_[1].first, T_conserved_variable_)
               .ViewComponent("cell", false);

  // in the case of water loss, use pressure.  in the case of water gain, use flux.
  for (int c = 0; c != p_star.MyLength(); ++c) {
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
  changedEvaluatorPrimary(p_lateral_flow_source_, tags_[1].second, *S_);
  changedEvaluatorPrimary(T_lateral_flow_source_, tags_[1].second, *S_);

  // mark p and WC as changed
  changedEvaluatorPrimary(p_primary_variable_, tags_[1].first, *S_);
  // changedEvaluatorPrimary(p_conserved_variable_, tags_[1].first, *S_);

  // copy to subsurface and mark as changed
  auto p_sub_owner = S_->GetRecord(p_sub_primary_variable_, tags_[1].first).owner();
  CopySurfaceToSubsurface(
    S_->Get<CompositeVector>(p_primary_variable_, tags_[1].first),
    S_->GetW<CompositeVector>(p_sub_primary_variable_, tags_[1].first, p_sub_owner));
  changedEvaluatorPrimary(p_sub_primary_variable_, tags_[1].first, *S_);

  // mark T and E as changed
  changedEvaluatorPrimary(T_primary_variable_, tags_[1].first, *S_);
  // changedEvaluatorPrimary(T_conserved_variable_, tags_[1].first, *S_);

  // copy to subsurface and mark as changed
  auto T_sub_owner = S_->GetRecord(T_sub_primary_variable_, tags_[1].first).owner();
  CopySurfaceToSubsurface(
    S_->Get<CompositeVector>(T_primary_variable_, tags_[1].first),
    S_->GetW<CompositeVector>(T_sub_primary_variable_, tags_[1].first, T_sub_owner));
  changedEvaluatorPrimary(T_sub_primary_variable_, tags_[1].first, *S_);
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system, using the pressure coupling schem
//
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_DomainSet_Pressure_()
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // copy p primary variables into star primary variable
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, tags_[0].second)
                          ->ViewComponent("cell", false);
  const auto& WC_star = *S_->GetPtr<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                           ->ViewComponent("cell", false);
  const auto& T_star = *S_->GetPtr<CompositeVector>(T_primary_variable_star_, tags_[0].second)
                          ->ViewComponent("cell", false);
  const auto& E_star = *S_->GetPtr<CompositeVector>(T_conserved_variable_star_, tags_[0].second)
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

      Key WC_key = Keys::getKey(*ds_iter, p_conserved_variable_suffix_);
      auto WC_owner = S_->GetRecord(WC_key, ds_tag_current).owner();
      auto& WC =
        *S_->GetW<CompositeVector>(WC_key, ds_tag_current, WC_owner).ViewComponent("cell", false);
      AMANZI_ASSERT(WC.MyLength() == 1);
      WC[0][0] = WC_star[0][c];

      changedEvaluatorPrimary(p_key, ds_tag_current, *S_);
      // changedEvaluatorPrimary(WC_key, ds_tag_current, *S_);

      Key p_sub_key = Keys::getKey(
        domain_sub_, Keys::getDomainSetIndex(*ds_iter), p_sub_primary_variable_suffix_);
      auto p_sub_owner = S_->GetRecord(p_sub_key, ds_tag_current).owner();
      CopySurfaceToSubsurface(S_->Get<CompositeVector>(p_key, ds_tag_current),
                              S_->GetW<CompositeVector>(p_sub_key, ds_tag_current, p_sub_owner));
    }

    Key T_key = Keys::getKey(*ds_iter, T_primary_variable_suffix_);
    auto T_owner = S_->GetRecord(T_key, ds_tag_current).owner();
    auto& T =
      *S_->GetW<CompositeVector>(T_key, ds_tag_current, T_owner).ViewComponent("cell", false);
    AMANZI_ASSERT(T.MyLength() == 1);
    T[0][0] = T_star[0][c];

    Key E_key = Keys::getKey(*ds_iter, T_conserved_variable_suffix_);
    auto E_owner = S_->GetRecord(E_key, ds_tag_current).owner();
    auto& E =
      *S_->GetW<CompositeVector>(E_key, ds_tag_current, E_owner).ViewComponent("cell", false);
    AMANZI_ASSERT(E.MyLength() == 1);
    E[0][0] = E_star[0][c];

    changedEvaluatorPrimary(T_key, ds_tag_current, *S_);
    // changedEvaluatorPrimary(E_key, ds_tag_current, *S_);

    Key T_sub_key =
      Keys::getKey(domain_sub_, Keys::getDomainSetIndex(*ds_iter), T_sub_primary_variable_suffix_);
    auto T_sub_owner = S_->GetRecord(T_sub_key, ds_tag_current).owner();
    CopySurfaceToSubsurface(S_->Get<CompositeVector>(T_key, ds_tag_current),
                            S_->GetW<CompositeVector>(T_sub_key, ds_tag_current, T_sub_owner));
    ++ds_iter;
  }
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_DomainSet_Flux_()
{
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // grab the data, difference
  Epetra_MultiVector q_div(*S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                              .ViewComponent("cell", false));
  q_div.Update(
    1.0 / dt,
    *S_->Get<CompositeVector>(p_conserved_variable_star_, tag_next_).ViewComponent("cell", false),
    -1.0 / dt,
    *S_->Get<CompositeVector>(p_conserved_variable_star_, tag_current_)
       .ViewComponent("cell", false),
    0.);
  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[0].second).ViewComponent("cell", false),
    q_div,
    0.);

  // grab the data, difference
  Epetra_MultiVector qE_div(*S_->Get<CompositeVector>(T_conserved_variable_star_, tags_[0].second)
                               .ViewComponent("cell", false));
  qE_div.Update(
    1.0 / dt,
    *S_->Get<CompositeVector>(T_conserved_variable_star_, tag_next_).ViewComponent("cell", false),
    -1.0 / dt,
    *S_->Get<CompositeVector>(T_conserved_variable_star_, tag_current_)
       .ViewComponent("cell", false),
    0.);
  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[0].second).ViewComponent("cell", false),
    qE_div,
    0.);

  // copy into columns
  auto ds_iter = domain_set.begin();
  for (int c = 0; c != q_div.MyLength(); ++c) {
    Tag ds_tag_next = get_ds_tag_next_(*ds_iter);
    Key p_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
    (*S_->GetW<CompositeVector>(p_key, ds_tag_next, name_).ViewComponent("cell", false))[0][0] =
      q_div[0][c];
    changedEvaluatorPrimary(p_key, ds_tag_next, *S_);

    Key T_key = Keys::getKey(*ds_iter, T_lateral_flow_source_suffix_);
    (*S_->GetW<CompositeVector>(T_key, ds_tag_next, name_).ViewComponent("cell", false))[0][0] =
      qE_div[0][c];
    changedEvaluatorPrimary(T_key, ds_tag_next, *S_);
    ++ds_iter;
  }
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_DomainSet_Hybrid_()
{
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // grab the data, difference
  Epetra_MultiVector q_div(*S_->Get<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                              .ViewComponent("cell", false));
  q_div.Update(
    1.0 / dt,
    *S_->Get<CompositeVector>(p_conserved_variable_star_, tag_next_).ViewComponent("cell", false),
    -1.0 / dt,
    *S_->Get<CompositeVector>(p_conserved_variable_star_, tag_current_)
       .ViewComponent("cell", false),
    0.);
  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[0].second).ViewComponent("cell", false),
    q_div,
    0.);

  // grab the data, difference
  Epetra_MultiVector qE_div(*S_->Get<CompositeVector>(T_conserved_variable_star_, tags_[0].second)
                               .ViewComponent("cell", false));
  qE_div.Update(
    1.0 / dt,
    *S_->Get<CompositeVector>(T_conserved_variable_star_, tag_next_).ViewComponent("cell", false),
    -1.0 / dt,
    *S_->Get<CompositeVector>(T_conserved_variable_star_, tag_current_)
       .ViewComponent("cell", false),
    0.);
  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(
    1.0,
    *S_->Get<CompositeVector>(cv_key_, tags_[0].second).ViewComponent("cell", false),
    qE_div,
    0.);

  // copy p primary variables into star primary variable
  const auto& p_star = *S_->GetPtr<CompositeVector>(p_primary_variable_star_, tags_[0].second)
                          ->ViewComponent("cell", false);
  const auto& WC_star = *S_->GetPtr<CompositeVector>(p_conserved_variable_star_, tags_[0].second)
                           ->ViewComponent("cell", false);
  const auto& T_star = *S_->GetPtr<CompositeVector>(T_primary_variable_star_, tags_[0].second)
                          ->ViewComponent("cell", false);
  const auto& E_star = *S_->GetPtr<CompositeVector>(T_conserved_variable_star_, tags_[0].second)
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

      Key WC_key = Keys::getKey(*ds_iter, p_conserved_variable_suffix_);
      auto WC_owner = S_->GetRecord(WC_key, ds_tag_current).owner();
      auto& WC =
        *S_->GetW<CompositeVector>(WC_key, ds_tag_current, WC_owner).ViewComponent("cell", false);
      AMANZI_ASSERT(WC.MyLength() == 1);
      WC[0][0] = WC_star[0][c];

      changedEvaluatorPrimary(p_key, ds_tag_current, *S_);
      // changedEvaluatorPrimary(WC_key, ds_tag_current, *S_);

      Key p_sub_key = Keys::getKey(
        domain_sub_, Keys::getDomainSetIndex(*ds_iter), p_sub_primary_variable_suffix_);
      auto p_sub_owner = S_->GetRecord(p_sub_key, ds_tag_current).owner();
      CopySurfaceToSubsurface(S_->Get<CompositeVector>(p_key, ds_tag_current),
                              S_->GetW<CompositeVector>(p_sub_key, ds_tag_current, p_sub_owner));

      Key T_key = Keys::getKey(*ds_iter, T_primary_variable_suffix_);
      auto T_owner = S_->GetRecord(T_key, ds_tag_current).owner();
      auto& T =
        *S_->GetW<CompositeVector>(T_key, ds_tag_current, T_owner).ViewComponent("cell", false);
      AMANZI_ASSERT(T.MyLength() == 1);
      T[0][0] = T_star[0][c];

      Key E_key = Keys::getKey(*ds_iter, T_conserved_variable_suffix_);
      auto E_owner = S_->GetRecord(E_key, ds_tag_current).owner();
      auto& E =
        *S_->GetW<CompositeVector>(E_key, ds_tag_current, E_owner).ViewComponent("cell", false);
      AMANZI_ASSERT(E.MyLength() == 1);
      E[0][0] = E_star[0][c];

      changedEvaluatorPrimary(T_key, ds_tag_current, *S_);
      // changedEvaluatorPrimary(E_key, ds_tag_current, *S_);

      Key T_sub_key = Keys::getKey(
        domain_sub_, Keys::getDomainSetIndex(*ds_iter), T_sub_primary_variable_suffix_);
      auto T_sub_owner = S_->GetRecord(T_sub_key, ds_tag_current).owner();
      CopySurfaceToSubsurface(S_->Get<CompositeVector>(T_key, ds_tag_current),
                              S_->GetW<CompositeVector>(T_sub_key, ds_tag_current, T_sub_owner));

      // set the lateral flux to 0
      Key p_lf_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
      (*S_->GetW<CompositeVector>(p_lf_key, ds_tag_next, name_)
          .ViewComponent("cell", false))[0][0] = 0.;
      changedEvaluatorPrimary(p_lf_key, ds_tag_next, *S_);

      Key T_lf_key = Keys::getKey(*ds_iter, T_lateral_flow_source_suffix_);
      (*S_->GetW<CompositeVector>(T_lf_key, ds_tag_next, name_)
          .ViewComponent("cell", false))[0][0] = 0.;
      changedEvaluatorPrimary(T_lf_key, ds_tag_next, *S_);

    } else {
      // use flux
      Key p_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
      (*S_->GetW<CompositeVector>(p_key, ds_tag_next, name_).ViewComponent("cell", false))[0][0] =
        q_div[0][c];
      changedEvaluatorPrimary(p_key, ds_tag_next, *S_);

      Key T_key = Keys::getKey(*ds_iter, T_lateral_flow_source_suffix_);
      (*S_->GetW<CompositeVector>(T_key, ds_tag_next, name_).ViewComponent("cell", false))[0][0] =
        qE_div[0][c];
      changedEvaluatorPrimary(T_key, ds_tag_next, *S_);
    }
    ++ds_iter;
  }
}


} // namespace Amanzi
