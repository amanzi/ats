/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"
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
  // collect domain names
  domain_set_ = Keys::readDomain(*plist_); // e.g. surface or surface_column:*
  domain_star_ = Keys::readDomain(*plist_, "star"); // e.g. surface_star

  // determine whether we are coupling subdomains or coupling 3D domains
  is_domain_set_ = S->HasDomainSet(domain_set_);
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

    p_conserved_variable_star_ = Keys::readKey(*plist_, domain_star_, "water conserved quantity star", "water_content");
    T_conserved_variable_star_ = Keys::readKey(*plist_, domain_star_, "energy conserved quantity star", "energy");

    cv_key_ = Keys::readKey(*plist_, domain_star_, "cell volume", "cell_volume");

    // set up for a primary variable field evaluator for the flux
    auto& p_sublist = S->FEList().sublist(p_lateral_flow_source_);
    p_sublist.set("field evaluator type", "primary variable");
    auto& T_sublist = S->FEList().sublist(T_lateral_flow_source_);
    T_sublist.set("field evaluator type", "primary variable");
  }

  // init sub-pks
  init_(S);

  // check whether we are subcycling
  subcycled_ = plist_->template get<bool>("subcycle subdomains", false);
  if (subcycled_) {
    subcycled_target_dt_ = plist_->template get<double>("subcycling target time step [s]");
    subcycled_min_dt_ = plist_->template get<double>("minimum subcycled time step [s]", 1.e-4);
  }
};


void MPCPermafrostSplitFlux::Initialize(const Teuchos::Ptr<State>& S)
{
  sub_pks_[1]->Initialize(S);
  CopyPrimaryToStar_(S, S);
  sub_pks_[0]->Initialize(S);

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
  S->GetField(p_primary_variable_star_, S->GetField(p_primary_variable_star_)->owner())->set_initialized();
  S->GetField(T_primary_variable_star_, S->GetField(T_primary_variable_star_)->owner())->set_initialized();

  auto& pstar = *S->GetFieldData(p_primary_variable_star_)->ViewComponent("cell", false);
  std::cout << "INITIALIZE: pstar = " << pstar[0][0] << std::endl;

  // set the fluxes as initialized -- they will get set by later calls to
  // CopyStarToPrimary, so no need for values, but do need to toggle the flag.
  if (coupling_ != "pressure") {
    if (is_domain_set_) {
      auto domain_set = S->GetDomainSet(domain_set_);
      for (const auto& domain : *domain_set) {
        auto pkey = Keys::getKey(domain, p_lateral_flow_source_suffix_);
        S->GetField(pkey, pkey)->set_initialized();

        auto Tkey = Keys::getKey(domain, T_lateral_flow_source_suffix_);
        S->GetField(Tkey, Tkey)->set_initialized();
      }
    } else {
      S->GetField(p_lateral_flow_source_, p_lateral_flow_source_)->set_initialized();
      S->GetField(T_lateral_flow_source_, T_lateral_flow_source_)->set_initialized();
    }
  }
}


void MPCPermafrostSplitFlux::Setup(const Teuchos::Ptr<State>& S)
{
  MPC<PK>::Setup(S);

  if (coupling_ != "pressure") {
    if (is_domain_set_) {
      auto domain_set = S->GetDomainSet(domain_set_);
      for (const auto& domain : *domain_set) {
        auto pkey = Keys::getKey(domain, p_lateral_flow_source_suffix_);
        S->RequireField(pkey, pkey)
          ->SetMesh(S->GetMesh(domain))
          ->SetComponent("cell", AmanziMesh::CELL, 1);
        S->RequireFieldEvaluator(pkey);

        auto Tkey = Keys::getKey(domain, T_lateral_flow_source_suffix_);
        S->RequireField(Tkey, Tkey)
          ->SetMesh(S->GetMesh(domain))
          ->SetComponent("cell", AmanziMesh::CELL, 1);
        S->RequireFieldEvaluator(Tkey);
      }
    } else {
      S->RequireField(p_lateral_flow_source_, p_lateral_flow_source_)
        ->SetMesh(S->GetMesh(Keys::getDomain(p_lateral_flow_source_)))
        ->SetComponent("cell", AmanziMesh::CELL, 1);
      S->RequireFieldEvaluator(p_lateral_flow_source_);

      S->RequireField(T_lateral_flow_source_, T_lateral_flow_source_)
        ->SetMesh(S->GetMesh(Keys::getDomain(T_lateral_flow_source_)))
        ->SetComponent("cell", AmanziMesh::CELL, 1);
      S->RequireFieldEvaluator(T_lateral_flow_source_);
    }

    // also need conserved quantities
    S->RequireField(p_conserved_variable_star_)
      ->SetMesh(S->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(p_conserved_variable_star_);

    S->RequireField(T_conserved_variable_star_)
      ->SetMesh(S->GetMesh(domain_star_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(T_conserved_variable_star_);
  }
}

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCPermafrostSplitFlux::get_dt()
{
  double dt = 1.0e99;
  if (subcycled_) {
    dt = subcycled_target_dt_;
  } else {
    for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
         pk != sub_pks_.end(); ++pk) {
      dt = std::min<double>(dt, (*pk)->get_dt());
    }
  }
  return dt;
};

// -----------------------------------------------------------------------------
// Set timestep for sub PKs
// -----------------------------------------------------------------------------
void MPCPermafrostSplitFlux::set_dt( double dt)
{
  if (subcycled_) {
    cycle_dt_ = dt;
  } else {
    for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
         pk != sub_pks_.end(); ++pk) {
      (*pk)->set_dt(dt);
    }
  }
};

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool MPCPermafrostSplitFlux::AdvanceStep(double t_old, double t_new, bool reinit)
{
  if (subcycled_) return AdvanceStep_Subcycled_(t_old, t_new, reinit);
  else return AdvanceStep_Standard_(t_old, t_new, reinit);
}


bool MPCPermafrostSplitFlux::AdvanceStep_Standard_(double t_old, double t_new, bool reinit)
{
  // Advance the star system
  bool fail = false;
  AMANZI_ASSERT(sub_pks_.size() == 2);
  fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;

  // Copy star's new value into primary's old value
  CopyStarToPrimary_(t_new - t_old);

  // Now advance the primary
  fail = sub_pks_[1]->AdvanceStep(t_old, t_new, reinit);
  return fail;
};


bool MPCPermafrostSplitFlux::AdvanceStep_Subcycled_(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  bool fail = false;
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Beginning subcycled timestepping." << std::endl;

  int my_pid = solution_->Comm()->MyPID();

  // note, both can fail, in the case of 3D subdomain.  In the case of columns,
  // the columns should be subcycled independently and therefore never fail,
  // but this won't break in that case.
  AMANZI_ASSERT(sub_pks_.size() == 2);
  for (int i=0; i!=sub_pks_.size(); ++i) {
    if (vo_->os_OK(Teuchos::VERB_EXTREME))
      *vo_->os() << "Beginning subcyling on pk \"" << sub_pks_[i]->name() << "\"" << std::endl;

    bool done = false;
    double t_inner = t_old;
    S_inter_->set_time(t_old);
    while (!done) {
      double dt_inner = std::min(sub_pks_[i]->get_dt(), t_new - t_inner);
      *S_next_->GetScalarData("dt", "coordinator") = dt_inner;
      S_next_->set_time(t_inner + dt_inner);
      bool fail_inner = sub_pks_[i]->AdvanceStep(t_inner, t_inner+dt_inner, false);
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  step failed? " << fail_inner << std::endl;
      bool valid_inner = sub_pks_[i]->ValidStep();
      if (vo_->os_OK(Teuchos::VERB_EXTREME))
        *vo_->os() << "  step valid? " << valid_inner << std::endl;

      if (fail_inner || !valid_inner) {
        // FIXME: figure out a way to get child domains, rather than guess at these! --etc
        if (i == 0) { // star system
          S_next_->AssignDomain(*S_inter_, domain_star_);
        } else {
          // note, this is always in the 3D case 
          S_next_->AssignDomain(*S_inter_, domain_sub_);
          S_next_->AssignDomain(*S_inter_, domain_);
          S_next_->AssignDomain(*S_inter_, domain_snow_);
        }

        dt_inner = sub_pks_[i]->get_dt();
        S_next_->set_time(S_inter_->time());
        S_next_->set_cycle(S_inter_->cycle());

        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  failed, new timestep is " << dt_inner << std::endl;

      } else {
        sub_pks_[i]->CommitStep(t_inner, t_inner + dt_inner, S_next_);
        t_inner += dt_inner;
        if (std::abs(t_new - t_inner) < 1.e-10) done = true;

        // FIXME: figure out a way to get child domains, rather than guess at these! --etc
        if (i == 0) { // star system
          S_inter_->AssignDomain(*S_inter_, domain_star_);
        } else {
          S_inter_->AssignDomain(*S_next_, domain_sub_);
          S_inter_->AssignDomain(*S_next_, domain_);
          S_inter_->AssignDomain(*S_next_, domain_snow_);
        }

        S_inter_->set_time(S_next_->time());
        S_inter_->set_cycle(S_next_->cycle());
        dt_inner = sub_pks_[i]->get_dt();
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          *vo_->os() << "  success, new timestep is " << dt_inner << std::endl;
      }

      if (dt_inner < subcycled_min_dt_) {
        Errors::Message msg;
        msg << "SubPK " << sub_pks_[i]->name() << " on PID " << my_pid << " crashing timestep in subcycling: dt = " << dt_inner;
        Exceptions::amanzi_throw(msg);
      }
    }

    // do the in-between things -- this part is not generic, the above is
    // generic enough for MPC!
    if (i == 0) CopyStarToPrimary_(t_new - t_old);
  }
  S_inter_->set_time(t_old);
  return false;
};


bool
MPCPermafrostSplitFlux::ValidStep() {
  if (subcycled_) return true; // this was already checked in advance
  else return MPC<PK>::ValidStep();
}


void MPCPermafrostSplitFlux::CommitStep(double t_old, double t_new,
        const Teuchos::RCP<State>& S)
{
  auto& pstar = *S->GetFieldData(p_primary_variable_star_)->ViewComponent("cell", false);

  // NOTE: in AJC code, these were flipped.  I believe this is correct, but
  // might be worth checking to see which works better.  Note it should result
  // in larget timestep sizes to get this right, but the physics shouldn't
  // change if this is wrong.  In AJC code, the below comment was still
  // present...
  //
  // Also, if this _did_ need to be flipped, than the call to
  // sub_pks_[i]->CommitStep() in the Advance_Subcycled would be incorrect, and
  // we would have to unravel that loop and not call Commit on the star system.
  if (!subcycled_) {
    // Commit before copy to ensure record for extrapolation in star system uses
    // its own solutions
    MPC<PK>::CommitStep(t_old, t_new, S);
  }

  // Copy the primary into the star to advance
  CopyPrimaryToStar_(S.ptr(), S.ptr());

  auto& pstar2 = *S->GetFieldData(p_primary_variable_star_)->ViewComponent("cell", false);
}


void
MPCPermafrostSplitFlux::CopyPrimaryToStar_(const Teuchos::Ptr<const State>& S,
                                    const Teuchos::Ptr<State>& S_star) {
  if (is_domain_set_) {
    CopyPrimaryToStar_DomainSet_(S, S_star);
  } else {
    CopyPrimaryToStar_Standard_(S, S_star);
  }
}

void
MPCPermafrostSplitFlux::CopyStarToPrimary_(double dt) {
  if (is_domain_set_) {
    if (coupling_ == "pressure") CopyStarToPrimary_DomainSet_Pressure_(dt);
    else if (coupling_ == "flux") CopyStarToPrimary_DomainSet_Flux_(dt);
    else if (coupling_ == "hybrid") CopyStarToPrimary_DomainSet_Hybrid_(dt);
    else AMANZI_ASSERT(false);
  } else {
    if (coupling_ == "pressure") CopyStarToPrimary_Standard_Pressure_(dt);
    else if (coupling_ == "flux") CopyStarToPrimary_Standard_Flux_(dt);
    else if (coupling_ == "hybrid") CopyStarToPrimary_Standard_Hybrid_(dt);
    else AMANZI_ASSERT(false);
  }
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system assuming a 3D subsurface domain
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyPrimaryToStar_Standard_(const Teuchos::Ptr<const State>& S,
                                    const Teuchos::Ptr<State>& S_star)
{
  // copy p primary variable
  auto& p_star = *S_star->GetFieldData(p_primary_variable_star_, S_star->GetField(p_primary_variable_star_)->owner())
                  ->ViewComponent("cell",false);
  const auto& p = *S->GetFieldData(p_primary_variable_)->ViewComponent("cell",false);
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p[0][c] <= 101325.0) {
      p_star[0][c] = 101325.;
    } else {
      p_star[0][c] = p[0][c];
    }
  }

  auto peval = S_star->GetFieldEvaluator(p_primary_variable_star_);
  auto peval_pvfe = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(peval);
  peval_pvfe->SetFieldAsChanged(S_star.ptr());

  // copy T primary variable
  auto& T_star = *S_star->GetFieldData(T_primary_variable_star_, S_star->GetField(T_primary_variable_star_)->owner())
                  ->ViewComponent("cell",false);
  const auto& T = *S->GetFieldData(T_primary_variable_)->ViewComponent("cell",false);
  T_star = T;

  auto Teval = S_star->GetFieldEvaluator(T_primary_variable_star_);
  auto Teval_pvfe = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(Teval);
  Teval_pvfe->SetFieldAsChanged(S_star.ptr());
}

// -----------------------------------------------------------------------------
// Copy the primary variable to the star system assuming a DomainSet domain
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyPrimaryToStar_DomainSet_(const Teuchos::Ptr<const State>& S,
                                    const Teuchos::Ptr<State>& S_star)
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // copy p primary variables into star primary variable
  auto& p_star = *S_star->GetFieldData(p_primary_variable_star_, S_star->GetField(p_primary_variable_star_)->owner())
                  ->ViewComponent("cell",false);
  auto ds_iter = domain_set.begin();
  for (int c=0; c!=p_star.MyLength(); ++c) {
    Key pkey = Keys::getKey(*ds_iter, p_primary_variable_suffix_);
    const auto& p = *S->GetFieldData(pkey)->ViewComponent("cell",false);
    AMANZI_ASSERT(p.MyLength() == 1);
    if (p[0][0] <= 101325.0) {
      p_star[0][c] = 101325.;
    } else {
      p_star[0][c] = p[0][0];
    }
    ++ds_iter;
  }
  auto peval = S_star->GetFieldEvaluator(p_primary_variable_star_);
  auto peval_pvfe = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(peval);
  peval_pvfe->SetFieldAsChanged(S_star.ptr());

  // copy T primary variable
  auto& T_star = *S_star->GetFieldData(T_primary_variable_star_, S_star->GetField(T_primary_variable_star_)->owner())
                  ->ViewComponent("cell",false);
  ds_iter = domain_set.begin();
  for (int c=0; c!=T_star.MyLength(); ++c) {
    Key Tkey = Keys::getKey(*ds_iter, T_primary_variable_suffix_);
    const auto& T = *S->GetFieldData(Tkey)->ViewComponent("cell",false);
    AMANZI_ASSERT(T.MyLength() == 1);
    T_star[0][c] = T[0][0];
    ++ds_iter;
  }
  auto Teval = S_star->GetFieldEvaluator(T_primary_variable_star_);
  auto Teval_pvfe = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(Teval);
  Teval_pvfe->SetFieldAsChanged(S_star.ptr());
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_Standard_Pressure_(double dt)
{
  // copy p primary variables from star primary variable
  const auto& p_star = *S_next_->GetFieldData(p_primary_variable_star_)
    ->ViewComponent("cell", false);
  auto passwd = S_inter_->GetField(p_primary_variable_)->owner();
  auto& p = *S_inter_->GetFieldData(p_primary_variable_, passwd)
    ->ViewComponent("cell", false);

  double p_atm_plus_eps = 101325. + 1.e-7;
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > p_atm_plus_eps) {
      p[0][c] = p_star[0][c];
    }
  }

  // mark p and subsurface p as changed
  {
    auto eval_p = S_inter_->GetFieldEvaluator(p_primary_variable_);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());
  }
  passwd = S_inter_->GetField(p_sub_primary_variable_)->owner();
  CopySurfaceToSubsurface(*S_inter_->GetFieldData(p_primary_variable_),
                          S_inter_->GetFieldData(p_sub_primary_variable_, passwd).ptr());
  {
    auto eval_p = S_inter_->GetFieldEvaluator(p_sub_primary_variable_);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());
  }


  // copy T primary variables from star primary variable
  const auto& T_star = *S_next_->GetFieldData(T_primary_variable_star_)
    ->ViewComponent("cell", false);
  passwd = S_inter_->GetField(T_primary_variable_)->owner();
  auto& T = *S_inter_->GetFieldData(T_primary_variable_, passwd)
    ->ViewComponent("cell", false);
  T = T_star;

  // mark T and subsurface T as changed
  {
    auto eval_p = S_inter_->GetFieldEvaluator(T_primary_variable_);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());
  }
  passwd = S_inter_->GetField(T_sub_primary_variable_)->owner();
  CopySurfaceToSubsurface(*S_inter_->GetFieldData(T_primary_variable_),
                          S_inter_->GetFieldData(T_sub_primary_variable_, passwd).ptr());
  {
    auto eval_p = S_inter_->GetFieldEvaluator(T_sub_primary_variable_);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());
  }
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_Standard_Flux_(double dt)
{
  // make sure we have the evaluator at the new state timestep
  if (p_eval_pvfe_ == Teuchos::null) {
    Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(p_lateral_flow_source_);
    p_eval_pvfe_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe);
    AMANZI_ASSERT(p_eval_pvfe_ != Teuchos::null);
  }
  if (T_eval_pvfe_ == Teuchos::null) {
    Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(T_lateral_flow_source_);
    T_eval_pvfe_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe);
    AMANZI_ASSERT(T_eval_pvfe_ != Teuchos::null);
  }

  // these updates should do nothing, but you never know
  S_inter_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);

  // grab the data, difference
  auto& q_div = *S_next_->GetFieldData(p_lateral_flow_source_, S_next_->GetField(p_lateral_flow_source_)->owner())
                ->ViewComponent("cell",false);
  q_div.Update(1.0/dt,
               *S_next_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), q_div, 0.);

  // mark the source evaluator as changed to ensure the total source gets updated.
  p_eval_pvfe_->SetFieldAsChanged(S_next_.ptr());

  // grab the data, difference
  auto& qE_div = *S_next_->GetFieldData(T_lateral_flow_source_, S_next_->GetField(T_lateral_flow_source_)->owner())
                ->ViewComponent("cell",false);
  qE_div.Update(1.0/dt,
               *S_next_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), qE_div, 0.);

  // mark the source evaluator as changed to ensure the total source gets updated.
  T_eval_pvfe_->SetFieldAsChanged(S_next_.ptr());
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_Standard_Hybrid_(double dt)
{
  // make sure we have the evaluator at the new state timestep
  if (p_eval_pvfe_ == Teuchos::null) {
    Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(p_lateral_flow_source_);
    p_eval_pvfe_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe);
    AMANZI_ASSERT(p_eval_pvfe_ != Teuchos::null);
  }
  if (T_eval_pvfe_ == Teuchos::null) {
    Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(T_lateral_flow_source_);
    T_eval_pvfe_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe);
    AMANZI_ASSERT(T_eval_pvfe_ != Teuchos::null);
  }

  // these updates should do nothing, but you never know
  S_inter_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);

  // grab the data, difference
  auto& q_div = *S_next_->GetFieldData(p_lateral_flow_source_, S_next_->GetField(p_lateral_flow_source_)->owner())
                ->ViewComponent("cell",false);
  q_div.Update(1.0/dt,
               *S_next_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), q_div, 0.);

  // grab the data, difference
  auto& qE_div = *S_next_->GetFieldData(T_lateral_flow_source_, S_next_->GetField(T_lateral_flow_source_)->owner())
                ->ViewComponent("cell",false);
  qE_div.Update(1.0/dt,
               *S_next_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), qE_div, 0.);

  // grab the pressure and temp from the star system as well
  const auto& p_star = *S_next_->GetFieldData(p_primary_variable_star_)
                       ->ViewComponent("cell",false);
  const auto& T_star = *S_next_->GetFieldData(T_primary_variable_star_)
                       ->ViewComponent("cell",false);

  // and from the surface system
  auto& p = *S_inter_->GetFieldData(p_primary_variable_,
          S_inter_->GetField(p_primary_variable_)->owner())
    ->ViewComponent("cell", false);
  auto& T = *S_inter_->GetFieldData(T_primary_variable_,
          S_inter_->GetField(T_primary_variable_)->owner())
    ->ViewComponent("cell", false);

  // in the case of water loss, use pressure.  in the case of water gain, use flux.
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > 101325. && q_div[0][c] < 0.) {
      // use the Dirichlet
      p[0][c] = p_star[0][c];
      T[0][c] = T_star[0][c];

      // set the lateral flux to 0
      q_div[0][c] = 0.;
      qE_div[0][c] = 0.;

    } // otherwise, use the flux, so nothing changes
  }

  // mark both the primary evals and flux evals as changed
  p_eval_pvfe_->SetFieldAsChanged(S_next_.ptr());
  T_eval_pvfe_->SetFieldAsChanged(S_next_.ptr());

  // mark p and subsurface p as changed
  {
    auto eval_p = S_inter_->GetFieldEvaluator(p_primary_variable_);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());
  }
  auto passwd = S_inter_->GetField(p_sub_primary_variable_)->owner();
  CopySurfaceToSubsurface(*S_inter_->GetFieldData(p_primary_variable_),
                          S_inter_->GetFieldData(p_sub_primary_variable_, passwd).ptr());
  {
    auto eval_p = S_inter_->GetFieldEvaluator(p_sub_primary_variable_);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());
  }
  // mark T and subsurface T as changed
  {
    auto eval_p = S_inter_->GetFieldEvaluator(T_primary_variable_);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());
  }
  passwd = S_inter_->GetField(T_sub_primary_variable_)->owner();
  CopySurfaceToSubsurface(*S_inter_->GetFieldData(T_primary_variable_),
                          S_inter_->GetFieldData(T_sub_primary_variable_, passwd).ptr());
  {
    auto eval_p = S_inter_->GetFieldEvaluator(T_sub_primary_variable_);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());
  }
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system, using the pressure coupling schem
//
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_DomainSet_Pressure_(double dt)
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // copy p primary variables into star primary variable
  const auto& p_star = *S_next_->GetFieldData(p_primary_variable_star_)
                       ->ViewComponent("cell",false);

  auto ds_iter = domain_set.begin();
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > 101325.0000001) {
      Key pkey = Keys::getKey(*ds_iter, p_primary_variable_suffix_);
      auto& p = *S_inter_->GetFieldData(pkey, S_inter_->GetField(pkey)->owner())->ViewComponent("cell",false);
      AMANZI_ASSERT(p.MyLength() == 1);
      p[0][0] = p_star[0][c];

      auto eval_p = S_inter_->GetFieldEvaluator(pkey);
      auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
      AMANZI_ASSERT(eval_pv.get());
      eval_pv->SetFieldAsChanged(S_inter_.ptr());

      Key psub_key = Keys::getKey(domain_sub_,
              Keys::getDomainSetIndex(*ds_iter), p_sub_primary_variable_suffix_);
      auto passwd = S_inter_->GetField(psub_key)->owner();
      CopySurfaceToSubsurface(*S_inter_->GetFieldData(pkey),
              S_inter_->GetFieldData(psub_key, passwd).ptr());
    }
    ++ds_iter;
  }

  // copy p primary variables into star primary variable
  ds_iter = domain_set.begin();
  const auto& T_star = *S_next_->GetFieldData(T_primary_variable_star_)
                       ->ViewComponent("cell",false);
  for (int c=0; c!=T_star.MyLength(); ++c) {
    Key Tkey = Keys::getKey(*ds_iter, T_primary_variable_suffix_);
    auto& T = *S_inter_->GetFieldData(Tkey, S_inter_->GetField(Tkey)->owner())->ViewComponent("cell",false);
    AMANZI_ASSERT(T.MyLength() == 1);
    T[0][0] = T_star[0][c];

    auto eval_p = S_inter_->GetFieldEvaluator(Tkey);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());

    Key Tsub_key = Keys::getKey(domain_sub_,
            Keys::getDomainSetIndex(*ds_iter), T_sub_primary_variable_suffix_);
    auto passwd = S_inter_->GetField(Tsub_key)->owner();
    CopySurfaceToSubsurface(*S_inter_->GetFieldData(Tkey),
                            S_inter_->GetFieldData(Tsub_key, passwd).ptr());
    ++ds_iter;
  }
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_DomainSet_Hybrid_(double dt)
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // make sure we have the evaluator at the new state timestep
  if (p_eval_pvfes_.size() == 0) {
    for (const auto& domain : domain_set) {
      Key p_lf_key = Keys::getKey(domain, p_lateral_flow_source_suffix_);
      Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(p_lf_key);
      p_eval_pvfes_.push_back(Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe));
      AMANZI_ASSERT(p_eval_pvfes_.back() != Teuchos::null);
    }
  }
  if (T_eval_pvfes_.size() == 0) {
    for (const auto& domain : domain_set) {
      Key T_lf_key = Keys::getKey(domain, T_lateral_flow_source_suffix_);
      Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(T_lf_key);
      T_eval_pvfes_.push_back(Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe));
      AMANZI_ASSERT(T_eval_pvfes_.back() != Teuchos::null);
    }
  }

  // these updates should do nothing, but you never know
  S_inter_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);

  // grab the data, difference
  Epetra_MultiVector q_div(*S_next_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false));
  q_div.Update(1.0/dt,
               *S_next_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), q_div, 0.);

  // grab the energy, difference
  Epetra_MultiVector qE_div(*S_next_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false));
  qE_div.Update(1.0/dt,
               *S_next_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               0.);
  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), qE_div, 0.);

  // grab the pressure from the star system as well
  const auto& p_star = *S_next_->GetFieldData(p_primary_variable_star_)
                       ->ViewComponent("cell",false);
  const auto& T_star = *S_next_->GetFieldData(T_primary_variable_star_)
                       ->ViewComponent("cell",false);

  // in the case of water loss, use pressure.  in the case of water gain, use flux.
  auto ds_iter = domain_set.begin();
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > 101325. && q_div[0][c] < 0.) {
      // use the Dirichlet
      Key pkey = Keys::getKey((*ds_iter), p_primary_variable_suffix_);
      auto& p = *S_inter_->GetFieldData(pkey, S_inter_->GetField(pkey)->owner())->ViewComponent("cell",false);
      AMANZI_ASSERT(p.MyLength() == 1);
      p[0][0] = p_star[0][c];

      Key Tkey = Keys::getKey(*ds_iter, T_primary_variable_suffix_);
      auto& T = *S_inter_->GetFieldData(Tkey, S_inter_->GetField(Tkey)->owner())->ViewComponent("cell",false);
      AMANZI_ASSERT(T.MyLength() == 1);
      T[0][0] = T_star[0][c];

      // tag the evaluators as changed
      auto eval_p = S_inter_->GetFieldEvaluator(pkey);
      auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
      AMANZI_ASSERT(eval_pv.get());
      eval_pv->SetFieldAsChanged(S_inter_.ptr());

      auto eval_t = S_inter_->GetFieldEvaluator(Tkey);
      auto eval_tv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_t);
      AMANZI_ASSERT(eval_tv.get());
      eval_tv->SetFieldAsChanged(S_inter_.ptr());

      // copy from surface to subsurface to ensure consistency
      Key psub_key = Keys::getKey(domain_sub_,
              Keys::getDomainSetIndex(*ds_iter), p_sub_primary_variable_suffix_);
      auto passwd = S_inter_->GetField(psub_key)->owner();
      CopySurfaceToSubsurface(*S_inter_->GetFieldData(pkey),
              S_inter_->GetFieldData(psub_key,passwd).ptr());
      Key Tsub_key = Keys::getKey(domain_sub_,
              Keys::getDomainSetIndex(*ds_iter), T_sub_primary_variable_suffix_);
      passwd = S_inter_->GetField(Tsub_key)->owner();
      CopySurfaceToSubsurface(*S_inter_->GetFieldData(Tkey),
              S_inter_->GetFieldData(Tsub_key, passwd).ptr());

      // set the lateral flux to 0
      Key p_lf_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
      (*S_next_->GetFieldData(p_lf_key, p_lf_key)->ViewComponent("cell",false))[0][0] = 0.;
      p_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());

      Key T_lf_key = Keys::getKey(*ds_iter, T_lateral_flow_source_suffix_);
      (*S_next_->GetFieldData(T_lf_key, T_lf_key)->ViewComponent("cell",false))[0][0] = 0.;
      T_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());

    } else {
      // use flux
      Key p_lf_key = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
      (*S_next_->GetFieldData(p_lf_key, p_lf_key)->ViewComponent("cell",false))[0][0] = q_div[0][c];
      p_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());

      Key T_lf_key = Keys::getKey(*ds_iter, T_lateral_flow_source_suffix_);
      (*S_next_->GetFieldData(T_lf_key, T_lf_key)->ViewComponent("cell",false))[0][0] = qE_div[0][c];
      T_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());
    }
    ++ds_iter;
  }
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFlux::CopyStarToPrimary_DomainSet_Flux_(double dt)
{
  const auto& domain_set = *S_->GetDomainSet(domain_set_);

  // make sure we have the evaluator at the new state timestep
  if (p_eval_pvfes_.size() == 0) {
    for (const auto& domain : domain_set) {
      Key pkey = Keys::getKey(domain, p_lateral_flow_source_suffix_);
      Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(pkey);
      p_eval_pvfes_.push_back(Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe));
      AMANZI_ASSERT(p_eval_pvfes_.back() != Teuchos::null);
    }
  }

  if (T_eval_pvfes_.size() == 0) {
    for (const auto& domain : domain_set) {
      Key pkey = Keys::getKey(domain, T_lateral_flow_source_suffix_);
      Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(pkey);
      T_eval_pvfes_.push_back(Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe));
      AMANZI_ASSERT(T_eval_pvfes_.back() != Teuchos::null);
    }
  }

  // these updates should do nothing, but you never know
  S_inter_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(p_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);
  S_inter_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_inter_.ptr(), name_);
  S_next_->GetFieldEvaluator(T_conserved_variable_star_)->HasFieldChanged(S_next_.ptr(), name_);

  // grab the data, difference
  Epetra_MultiVector q_div(*S_next_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false));
  q_div.Update(1.0/dt,
               *S_next_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(p_conserved_variable_star_)->ViewComponent("cell",false),
               0.);
  // scale by cell volume as this will get rescaled in the source calculation
  q_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), q_div, 0.);

  // copy into columns
  auto ds_iter = domain_set.begin();
  for (int c=0; c!=q_div.MyLength(); ++c) {
    Key pkey = Keys::getKey(*ds_iter, p_lateral_flow_source_suffix_);
    (*S_next_->GetFieldData(pkey, pkey)->ViewComponent("cell",false))[0][0] = q_div[0][c];
    p_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());
    ++ds_iter;
  }

  // grab the data, difference
  Epetra_MultiVector qE_div(*S_next_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false));
  qE_div.Update(1.0/dt,
               *S_next_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               -1.0/dt,
               *S_inter_->GetFieldData(T_conserved_variable_star_)->ViewComponent("cell",false),
               0.);

  // scale by cell volume as this will get rescaled in the source calculation
  qE_div.ReciprocalMultiply(1.0, *S_next_->GetFieldData(cv_key_)->ViewComponent("cell",false), qE_div, 0.);

  // copy into columns
  ds_iter = domain_set.begin();
  for (int c=0; c!=qE_div.MyLength(); ++c) {
    Key Tkey = Keys::getKey(*ds_iter, T_lateral_flow_source_suffix_);
    (*S_next_->GetFieldData(Tkey, Tkey)->ViewComponent("cell",false))[0][0] = qE_div[0][c];
    T_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());
    ++ds_iter;
  }
}


} // namespace


