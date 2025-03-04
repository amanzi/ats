/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Process kernel for coupling of Transport PK and Chemistry PK.
*/

#include "mpc_reactivetransport.hh"
#include "pk_helpers.hh"
#include "chem_pk_helpers.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCReactiveTransport::MPCReactiveTransport(Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                           const Teuchos::RCP<State>& S,
                                           const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_list, S, soln), WeakMPC(pk_tree, global_list, S, soln),
    chem_step_succeeded_(true),
    trans_step_succeeded_(true)
{
  alquimia_timer_ = Teuchos::TimeMonitor::getNewCounter("alquimia " + name());
}


void
MPCReactiveTransport::parseParameterList()
{
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string>>("PKs order");
  domain_ = getSubPKPlist_(1)->get<std::string>("domain name", "domain");

  // chemistry and transport share the same primary variable, we take it from transport
  tcc_key_ = Keys::readKey(*getSubPKPlist_(1), domain_, "primary variable key",
                           "total_component_concentration");
  getSubPKPlist_(0)->set<std::string>("primary variable key", tcc_key_);

  // both pks need access to the primary variable
  getSubPKPlist_(0)->set<std::string>("primary variable password", name_);
  getSubPKPlist_(1)->set<std::string>("primary variable password", name_);

  // tell chemistry to operator split
  getSubPKPlist_(0)->set("operator split", true);

  // Only one PK needs to set the IC -- we use chemistry to process geochemical
  // initial conditions -- but the "initial conditions" list must be present
  // for all PK_Physical PKs, so just touch it here to make sure it exists.
  getSubPKPlist_(1)->sublist("initial conditions");

  mol_dens_key_ = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");

  cast_sub_pks_();
  chemistry_pk_->parseParameterList();

#ifdef ALQUIMIA_ENABLED
  transport_pk_->setChemEngine(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_));
#endif

  // now transport can be parse with knowledge of names
  transport_pk_->parseParameterList();
}


void
MPCReactiveTransport::Setup()
{
  transport_pk_->Setup();
  chemistry_pk_->Setup();

  S_->Require<CompositeVector, CompositeVectorSpace>(mol_dens_key_, tag_next_)
    .SetMesh(S_->GetMesh(domain_))
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(mol_dens_key_, tag_next_);
}

void
MPCReactiveTransport::cast_sub_pks_()
{
  // cast and call parse on chemistry
  chemistry_pk_ = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(sub_pks_[0]);
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);

  // now chem engine is set and we can hand it to transport
  transport_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(sub_pks_[1]);
  AMANZI_ASSERT(transport_pk_ != Teuchos::null);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void
MPCReactiveTransport::Initialize()
{
  // NOTE: this requires that Reactive-Transport is done last, or at least
  // after the density of water can be evaluated.  This could be problematic
  // for, e.g., salinity intrusion problems where water density is a function
  // of concentration itself, but should work for all other problems?
  Teuchos::RCP<Epetra_MultiVector> tcc =
    S_->GetW<CompositeVector>(tcc_key_, tag_next_, name_).ViewComponent("cell", true);
  S_->GetEvaluator(mol_dens_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> mol_dens =
    S_->Get<CompositeVector>(mol_dens_key_, tag_next_).ViewComponent("cell", true);

  convertConcentrationToAmanzi(*mol_dens, *tcc, *tcc);
  chemistry_pk_->Initialize();
  convertConcentrationToATS(*mol_dens, *tcc, *tcc);
  transport_pk_->Initialize();
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
MPCReactiveTransport::get_dt()
{
  double dt_tran = transport_pk_->get_dt();
  double dt_chem = chemistry_pk_->get_dt();

  // chemistry should manage its own timestep, and which point the default
  // MPC::get_dt() is fine.  --ETC
  if (trans_step_succeeded_ && !chem_step_succeeded_) dt_chem *= 0.5;
  return std::min(dt_tran, dt_chem);
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
MPCReactiveTransport::AdvanceStep(double t_old, double t_new, bool reinit)
{
  trans_step_succeeded_ = false;
  chem_step_succeeded_ = false;

  // First we do a transport step.
  bool fail = transport_pk_->AdvanceStep(t_old, t_new, reinit);
  if (fail) return fail;
  trans_step_succeeded_ = true;

  // Second, we do a chemistry step.
  int num_aq_componets = transport_pk_->get_num_aqueous_component();

  Teuchos::RCP<Epetra_MultiVector> tcc =
    S_->GetW<CompositeVector>(tcc_key_, tag_next_, name_).ViewComponent("cell", true);

  S_->GetEvaluator(mol_dens_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> mol_dens =
    S_->Get<CompositeVector>(mol_dens_key_, tag_next_).ViewComponent("cell", true);

  fail |=
    advanceChemistry(*chemistry_pk_, t_old, t_new, reinit, *mol_dens, *tcc, *alquimia_timer_);
  changedEvaluatorPrimary(tcc_key_, tag_next_, *S_);
  if (!fail) chem_step_succeeded_ = true;

  return fail;
};


} // namespace Amanzi
