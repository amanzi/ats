/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "primary_variable_field_evaluator.hh"
#include "mpc_surface_subsurface_helpers.hh"

#include "mpc_permafrost_split_flux_columns.hh"

#include "PK_Physical.hh"

namespace Amanzi {

MPCPermafrostSplitFluxColumns::MPCPermafrostSplitFluxColumns(Teuchos::ParameterList& FElist,
                 const Teuchos::RCP<Teuchos::ParameterList>& plist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution),
      MPC<PK>(FElist, plist, S, solution)
{
  // collect keys and names
  std::string domain_star = plist_->get<std::string>("star domain name");

  // generate a list of sub-pks, based upon the star system + a list of columns
  // -- get the original PKs order -- the column must be the last one
  auto subpks = plist_->get<Teuchos::Array<std::string> >("PKs order");
  std::string colname = subpks[subpks.size()-1];
  subpks.pop_back();

  // subcycle or not
  coupling_ = plist_->get<std::string>("coupling type", "pressure");
  if (coupling_ != "pressure" && coupling_ != "flux" && coupling_ != "hybrid") {
    Errors::Message msg("WeakMPCSemiCoupled: \"coupling type\" must be one of \"pressure\", \"flux\", or \"hybrid\".");
    Exceptions::amanzi_throw(msg);
  }

  // -- generate the column triple domain set
  KeyTriple col_triple;
  bool is_ds = Keys::splitDomainSet(colname, col_triple);
  if (!is_ds) {
    Errors::Message msg;
    msg << "WeakMPCSemiCoupled subpk: \"" << colname << "\" should be a domain-set PK of the form column_*-NAME";
    Exceptions::amanzi_throw(msg);
  }

  std::string domain_col = std::get<0>(col_triple);
  std::string domain_surf = "surface_"+domain_col;
  
  // -- add for the various columns based on GIDs of the surface system
  auto surf_mesh = S->GetMesh(domain_star);
  int ncols = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int i=0; i!=ncols; ++i) {
    int gid = surf_mesh->cell_map(false).GID(i);
    std::stringstream domain_name_stream;
    domain_name_stream << domain_col << "_" << gid;
    subpks.push_back(Keys::getKey(domain_name_stream.str(), std::get<2>(col_triple)));
    col_domains_.push_back(domain_name_stream.str());
  }

  // set up keys
  p_primary_variable_suffix_ = plist_->get<std::string>("pressure primary variable suffix", "pressure");
  T_primary_variable_suffix_ = plist_->get<std::string>("temperature primary variable suffix", "temperature");

  p_primary_variable_star_ = Keys::readKey(*plist_, domain_star, "pressure primary variable star", Keys::getVarName(p_primary_variable_suffix_));
  T_primary_variable_star_ = Keys::readKey(*plist_, domain_star, "temperature primary variable star", Keys::getVarName(T_primary_variable_suffix_));

  if (coupling_ != "pressure") {
    p_lateral_flow_source_suffix_ = plist_->get<std::string>("mass lateral flow source suffix", "mass_lateral_flow_source");
    T_lateral_flow_source_suffix_ = plist_->get<std::string>("energy lateral flow source suffix", "energy_lateral_flow_source");

    p_conserved_variable_star_ = Keys::readKey(*plist_, domain_star, "mass conserved quantity star", "water_content");
    T_conserved_variable_star_ = Keys::readKey(*plist_, domain_star, "energy conserved quantity star", "energy");

    cv_key_ = Keys::readKey(*plist_, domain_star, "cell volume", "cell_volume");
  
    // set up for a primary variable field evaluator for the flux
    Key primary_pkey = Keys::getKey(domain_surf + std::string("_*"), p_lateral_flow_source_suffix_);
    auto& p_sublist = S->FEList().sublist(primary_pkey);
    p_sublist.set("field evaluator type", "primary variable");
    Key primary_Tkey = Keys::getKey(domain_surf + std::string("_*"), T_lateral_flow_source_suffix_);
    auto& T_sublist = S->FEList().sublist(primary_Tkey);
    T_sublist.set("field evaluator type", "primary variable");
  }

  // init sub-pks
  plist_->set("PKs order", subpks);
  init_(S);
};


void MPCPermafrostSplitFluxColumns::Initialize(const Teuchos::Ptr<State>& S)
{
  // initialize the columns
  for (int i=1; i!=sub_pks_.size(); ++i) sub_pks_[i]->Initialize(S);

  // copy the columns to the star system, which initializes the star system consistently without an IC
  CopyPrimaryToStar(S, S);

  // initialize the star system (IC already set, but other things might happen)
  sub_pks_[0]->Initialize(S);
  S->GetField(p_primary_variable_star_, S->GetField(p_primary_variable_star_)->owner())->set_initialized();
  S->GetField(T_primary_variable_star_, S->GetField(T_primary_variable_star_)->owner())->set_initialized();
}


void MPCPermafrostSplitFluxColumns::Setup(const Teuchos::Ptr<State>& S)
{
  MPC<PK>::Setup(S);

  if (coupling_ != "pressure") {
    // require the coupling sources
    for (const auto& col_domain : col_domains_) {
      std::string col_surf_domain = "surface_" + col_domain;
      Key pkey = Keys::getKey(col_surf_domain, p_lateral_flow_source_suffix_);
      S->RequireField(pkey, pkey)
          ->SetMesh(S->GetMesh(col_surf_domain))
          ->SetComponent("cell", AmanziMesh::CELL, 1);
      S->RequireFieldEvaluator(pkey);

      Key Tkey = Keys::getKey(col_surf_domain, T_lateral_flow_source_suffix_);
      S->RequireField(Tkey, Tkey)
          ->SetMesh(S->GetMesh(col_surf_domain))
          ->SetComponent("cell", AmanziMesh::CELL, 1);
      S->RequireFieldEvaluator(Tkey);
    }
  }
}

// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCPermafrostSplitFluxColumns::get_dt()
{
  double dt_l = 1.e99;
  for (auto pk : sub_pks_) {
    dt_l = std::min(pk->get_dt(), dt_l);
  }    
  double dt_g;
  S_next_->GetMesh(Keys::getDomain(p_primary_variable_star_))->get_comm()->MinAll(&dt_l, &dt_g, 1);
  return dt_g;
};

// -----------------------------------------------------------------------------
// Set timestep for sub PKs 
// -----------------------------------------------------------------------------
void MPCPermafrostSplitFluxColumns::set_dt( double dt)
{
  AMANZI_ASSERT(false);
  sub_pks_[0]->set_dt(dt);
  for (MPC<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt);
  }
};

// -----------------------------------------------------------------------------
// Advance each sub-PK individually.
// -----------------------------------------------------------------------------
bool MPCPermafrostSplitFluxColumns::AdvanceStep(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  // Advance the star system 
  bool fail = false;
  fail = sub_pks_[0]->AdvanceStep(t_old, t_new, reinit);
  fail |= !sub_pks_[0]->ValidStep();
  if (fail) return fail;

  // Copy star's new value into primary's old value
  CopyStarToPrimary(t_new - t_old);

  // Now advance the primary
  for (int i=1; i!=sub_pks_.size(); ++i) {
    fail |= sub_pks_[i]->AdvanceStep(t_old, t_new, reinit);
    if (fail) break;
    fail |= !sub_pks_[i]->ValidStep();
    if (fail) break;
  }

  int fail_l(fail);
  int fail_g;
  S_next_->GetMesh(Keys::getDomain(p_primary_variable_star_))->get_comm()->MaxAll(&fail_l, &fail_g, 1);
  return fail_g > 0;
};


bool MPCPermafrostSplitFluxColumns::ValidStep() 
{
  return MPC<PK>::ValidStep();
}


void MPCPermafrostSplitFluxColumns::CommitStep(double t_old, double t_new,
        const Teuchos::RCP<State>& S)
{

  // Copy the primary into the star to advance
  CopyPrimaryToStar(S.ptr(), S.ptr());

  // commit before copy to ensure record for extrapolation in star system uses
  // its own solutions
  MPC<PK>::CommitStep(t_old, t_new, S);


}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFluxColumns::CopyPrimaryToStar(const Teuchos::Ptr<const State>& S,
                                    const Teuchos::Ptr<State>& S_star)
{
  // copy p primary variables into star primary variable
  auto& p_star = *S_star->GetFieldData(p_primary_variable_star_, S_star->GetField(p_primary_variable_star_)->owner())
                  ->ViewComponent("cell",false);
  for (int c=0; c!=p_star.MyLength(); ++c) {
    Key pkey = Keys::getKey("surface_"+col_domains_[c], p_primary_variable_suffix_);
    const auto& p = *S->GetFieldData(pkey)->ViewComponent("cell",false);
    AMANZI_ASSERT(p.MyLength() == 1);
    if (p[0][0] <= 101325.0) {
      p_star[0][c] = 101325.;
    } else {
      p_star[0][c] = p[0][0];
    }
  }

  auto peval = S_star->GetFieldEvaluator(p_primary_variable_star_);
  auto peval_pvfe = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(peval);
  peval_pvfe->SetFieldAsChanged(S_star.ptr());

  // copy T primary variable
  auto& T_star = *S_star->GetFieldData(T_primary_variable_star_, S_star->GetField(T_primary_variable_star_)->owner())
                  ->ViewComponent("cell",false);
  for (int c=0; c!=p_star.MyLength(); ++c) {
    Key Tkey = Keys::getKey("surface_"+col_domains_[c], T_primary_variable_suffix_);
    const auto& T = *S->GetFieldData(Tkey)->ViewComponent("cell",false);
    AMANZI_ASSERT(T.MyLength() == 1);
    T_star[0][c] = T[0][0];
  }

  auto Teval = S_star->GetFieldEvaluator(T_primary_variable_star_);
  auto Teval_pvfe = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(Teval);
  Teval_pvfe->SetFieldAsChanged(S_star.ptr());
}


void 
MPCPermafrostSplitFluxColumns::CopyStarToPrimary(double dt) {
  if (coupling_ == "pressure") CopyStarToPrimaryPressure_(dt);
  else if (coupling_ == "flux") CopyStarToPrimaryFlux_(dt);
  else if (coupling_ == "hybrid") CopyStarToPrimaryHybrid_(dt);
  else AMANZI_ASSERT(false);
}

// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFluxColumns::CopyStarToPrimaryPressure_(double dt)
{
  // copy p primary variables into star primary variable
  const auto& p_star = *S_next_->GetFieldData(p_primary_variable_star_)
                       ->ViewComponent("cell",false);
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > 101325.0000001) {
      Key pkey = Keys::getKey("surface_"+col_domains_[c], p_primary_variable_suffix_);
      auto& p = *S_inter_->GetFieldData(pkey, S_inter_->GetField(pkey)->owner())->ViewComponent("cell",false);
      AMANZI_ASSERT(p.MyLength() == 1);
      p[0][0] = p_star[0][c];

      auto eval_p = S_inter_->GetFieldEvaluator(pkey);
      auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
      AMANZI_ASSERT(eval_pv.get());
      eval_pv->SetFieldAsChanged(S_inter_.ptr());

      CopySurfaceToSubsurface(*S_inter_->GetFieldData(pkey),
              S_inter_->GetFieldData(Keys::getKey(col_domains_[c], p_primary_variable_suffix_),
                      S_inter_->GetField(Keys::getKey(col_domains_[c], p_primary_variable_suffix_))->owner()).ptr());
    }
  }

  // copy p primary variables into star primary variable
  const auto& T_star = *S_next_->GetFieldData(T_primary_variable_star_)
                       ->ViewComponent("cell",false);
  for (int c=0; c!=T_star.MyLength(); ++c) {
    Key Tkey = Keys::getKey("surface_"+col_domains_[c], T_primary_variable_suffix_);
    auto& T = *S_inter_->GetFieldData(Tkey, S_inter_->GetField(Tkey)->owner())->ViewComponent("cell",false);
    AMANZI_ASSERT(T.MyLength() == 1);
    T[0][0] = T_star[0][c];

    auto eval_p = S_inter_->GetFieldEvaluator(Tkey);
    auto eval_pv = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(eval_p);
    AMANZI_ASSERT(eval_pv.get());
    eval_pv->SetFieldAsChanged(S_inter_.ptr());
    
    CopySurfaceToSubsurface(*S_inter_->GetFieldData(Tkey),
                            S_inter_->GetFieldData(Keys::getKey(col_domains_[c], T_primary_variable_suffix_),
                                    S_inter_->GetField(Keys::getKey(col_domains_[c], T_primary_variable_suffix_))->owner()).ptr());
  }
}


// -----------------------------------------------------------------------------
// Copy the primary variable to the star system
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFluxColumns::CopyStarToPrimaryHybrid_(double dt)
{
  // make sure we have the evaluator at the new state timestep
  if (p_eval_pvfes_.size() == 0) {
    for (const auto& col_domain : col_domains_) {
      Key p_lf_key = Keys::getKey("surface_"+col_domain, p_lateral_flow_source_suffix_);
      Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(p_lf_key);
      p_eval_pvfes_.push_back(Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe));
      AMANZI_ASSERT(p_eval_pvfes_.back() != Teuchos::null);
    }
  }

  if (T_eval_pvfes_.size() == 0) {
    for (const auto& col_domain : col_domains_) {
      Key T_lf_key = Keys::getKey("surface_"+col_domain, T_lateral_flow_source_suffix_);
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
  for (int c=0; c!=p_star.MyLength(); ++c) {
    if (p_star[0][c] > 101325. && q_div[0][c] < 0.) {
      // use the Dirichlet
      Key pkey = Keys::getKey("surface_"+col_domains_[c], p_primary_variable_suffix_);
      auto& p = *S_inter_->GetFieldData(pkey, S_inter_->GetField(pkey)->owner())->ViewComponent("cell",false);
      AMANZI_ASSERT(p.MyLength() == 1);
      p[0][0] = p_star[0][c];

      Key Tkey = Keys::getKey("surface_"+col_domains_[c], T_primary_variable_suffix_);
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
      CopySurfaceToSubsurface(*S_inter_->GetFieldData(pkey),
                              S_inter_->GetFieldData(Keys::getKey(col_domains_[c], p_primary_variable_suffix_),
                                                     S_inter_->GetField(Keys::getKey(col_domains_[c], p_primary_variable_suffix_))->owner()).ptr());
      CopySurfaceToSubsurface(*S_inter_->GetFieldData(Tkey),
                              S_inter_->GetFieldData(Keys::getKey(col_domains_[c], T_primary_variable_suffix_),
                                                     S_inter_->GetField(Keys::getKey(col_domains_[c], T_primary_variable_suffix_))->owner()).ptr());

      // set the lateral flux to 0
      Key p_lf_key = Keys::getKey("surface_"+col_domains_[c], p_lateral_flow_source_suffix_);
      (*S_next_->GetFieldData(p_lf_key, p_lf_key)->ViewComponent("cell",false))[0][0] = 0.;
      p_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());

      Key T_lf_key = Keys::getKey("surface_"+col_domains_[c], T_lateral_flow_source_suffix_);
      (*S_next_->GetFieldData(T_lf_key, T_lf_key)->ViewComponent("cell",false))[0][0] = 0.;
      T_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());

    } else { 
      // use flux
      Key p_lf_key = Keys::getKey("surface_"+col_domains_[c], p_lateral_flow_source_suffix_);
      (*S_next_->GetFieldData(p_lf_key, p_lf_key)->ViewComponent("cell",false))[0][0] = q_div[0][c];
      p_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());

      Key T_lf_key = Keys::getKey("surface_"+col_domains_[c], T_lateral_flow_source_suffix_);
      (*S_next_->GetFieldData(T_lf_key, T_lf_key)->ViewComponent("cell",false))[0][0] = qE_div[0][c];
      T_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());
    }
  }
}


// -----------------------------------------------------------------------------
// Copy the star time derivative to the source evaluator.
// -----------------------------------------------------------------------------
void
MPCPermafrostSplitFluxColumns::CopyStarToPrimaryFlux_(double dt)
{
  // make sure we have the evaluator at the new state timestep
  if (p_eval_pvfes_.size() == 0) {
    for (const auto& col_domain : col_domains_) {
      Key pkey = Keys::getKey("surface_"+col_domain, p_lateral_flow_source_suffix_);
      Teuchos::RCP<FieldEvaluator> fe = S_next_->GetFieldEvaluator(pkey);
      p_eval_pvfes_.push_back(Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fe));
      AMANZI_ASSERT(p_eval_pvfes_.back() != Teuchos::null);
    }
  }
  
  if (T_eval_pvfes_.size() == 0) {
    for (const auto& col_domain : col_domains_) {
      Key pkey = Keys::getKey("surface_"+col_domain, T_lateral_flow_source_suffix_);
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
  for (int c=0; c!=q_div.MyLength(); ++c) {
    Key pkey = Keys::getKey("surface_"+col_domains_[c], p_lateral_flow_source_suffix_);
    (*S_next_->GetFieldData(pkey, pkey)->ViewComponent("cell",false))[0][0] = q_div[0][c];
    p_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());
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
  for (int c=0; c!=qE_div.MyLength(); ++c) {
    Key Tkey = Keys::getKey("surface_"+col_domains_[c], T_lateral_flow_source_suffix_);
    (*S_next_->GetFieldData(Tkey, Tkey)->ViewComponent("cell",false))[0][0] = qE_div[0][c];
    T_eval_pvfes_[c]->SetFieldAsChanged(S_next_.ptr());
  }
}

// protected constructor of subpks
void MPCPermafrostSplitFluxColumns::init_(const Teuchos::RCP<State>& S)
{
  PKFactory pk_factory;
  Teuchos::Array<std::string> pk_order = plist_->get< Teuchos::Array<std::string> >("PKs order");
  int npks = pk_order.size();

  // create the coupling PK
  // -- create the solution vector
  auto pk_soln = Teuchos::rcp(new TreeVector());
  solution_->PushBack(pk_soln);

  // -- create the PK
  sub_pks_.push_back(pk_factory.CreatePK(pk_order[0], pk_tree_, global_list_, S, pk_soln));

  // create the column PKs
  for (int i=1; i!=npks; ++i) {
    // create the solution vector
    //auto pk_soln = Teuchos::rcp(new TreeVector(Epetra_MpiComm(MPI_COMM_SELF)));
    auto pk_soln = Teuchos::rcp(new TreeVector());
    solution_->PushBack(pk_soln);

    // create the PK
    std::string name_i = pk_order[i];
    sub_pks_.push_back(pk_factory.CreatePK(name_i, pk_tree_, global_list_, S, pk_soln));
  }
};

} // namespace
