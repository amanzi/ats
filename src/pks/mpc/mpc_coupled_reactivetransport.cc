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

#include "mpc_coupled_reactivetransport.hh"
#include "pk_helpers.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCCoupledReactiveTransport::MPCCoupledReactiveTransport(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& global_list,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_list, S, soln), WeakMPC(pk_tree, global_list, S, soln)
{
  chem_step_succeeded_ = true;

  alquimia_timer_ = Teuchos::TimeMonitor::getNewCounter("alquimia " + name());
  alquimia_surf_timer_ = Teuchos::TimeMonitor::getNewCounter("alquimia surface " + name());

  domain_ = Keys::readDomain(*plist_, "domain", "domain");
  domain_surf_ = Keys::readDomainHint(*plist_, domain_, "domain", "surface");

  tcc_key_ = Keys::readKey(
    *plist_, domain_, "total component concentration", "total_component_concentration");
  tcc_surf_key_ = Keys::readKey(*plist_,
                                domain_surf_,
                                "surface total component concentration",
                                "total_component_concentration");
  mol_dens_key_ = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");
  mol_dens_surf_key_ =
    Keys::readKey(*plist_, domain_surf_, "surface molar density liquid", "molar_density_liquid");
}


void
MPCCoupledReactiveTransport::Setup()
{
  cast_sub_pks_();

  // must Setup transport first to get alias for saturation, etc set up correctly
  coupled_transport_pk_->Setup();
  coupled_chemistry_pk_->Setup();

  S_->Require<CompositeVector, CompositeVectorSpace>(tcc_key_, tag_next_, "state")
    .SetMesh(S_->GetMesh(domain_))
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, chemistry_pk_->num_aqueous_components());
  S_->RequireEvaluator(tcc_key_, tag_next_);

  S_->Require<CompositeVector, CompositeVectorSpace>(tcc_surf_key_, tag_next_, "state")
    .SetMesh(S_->GetMesh(domain_surf_))
    ->SetGhosted()
    ->AddComponent(
      "cell", AmanziMesh::Entity_kind::CELL, chemistry_pk_surf_->num_aqueous_components());
  S_->RequireEvaluator(tcc_surf_key_, tag_next_);

  S_->Require<CompositeVector, CompositeVectorSpace>(mol_dens_key_, tag_next_)
    .SetMesh(S_->GetMesh(domain_))
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(mol_dens_key_, tag_next_);

  S_->Require<CompositeVector, CompositeVectorSpace>(mol_dens_surf_key_, tag_next_)
    .SetMesh(S_->GetMesh(domain_surf_))
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(mol_dens_surf_key_, tag_next_);
}


void
MPCCoupledReactiveTransport::cast_sub_pks_()
{
  coupled_transport_pk_ = Teuchos::rcp_dynamic_cast<MPCCoupledTransport>(sub_pks_[1]);
  AMANZI_ASSERT(coupled_transport_pk_ != Teuchos::null);

  coupled_chemistry_pk_ = Teuchos::rcp_dynamic_cast<WeakMPC>(sub_pks_[0]);
  AMANZI_ASSERT(coupled_chemistry_pk_ != Teuchos::null);

  transport_pk_ =
    Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(coupled_transport_pk_->get_subpk(0));
  AMANZI_ASSERT(transport_pk_ != Teuchos::null);
  transport_pk_surf_ =
    Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(coupled_transport_pk_->get_subpk(1));
  AMANZI_ASSERT(transport_pk_surf_ != Teuchos::null);

  chemistry_pk_ =
    Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(coupled_chemistry_pk_->get_subpk(0));
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);
  chemistry_pk_surf_ =
    Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(coupled_chemistry_pk_->get_subpk(1));
  AMANZI_ASSERT(chemistry_pk_surf_ != Teuchos::null);

  AMANZI_ASSERT(transport_pk_->domain() == chemistry_pk_->domain());
  AMANZI_ASSERT(transport_pk_surf_->domain() == chemistry_pk_surf_->domain());

// communicate chemistry engine to transport.
#ifdef ALQUIMIA_ENABLED
  transport_pk_->SetupAlquimia(
    Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_),
    chemistry_pk_->chem_engine());
  transport_pk_surf_->SetupAlquimia(
    Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_surf_),
    chemistry_pk_surf_->chem_engine());
#endif
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void
MPCCoupledReactiveTransport::Initialize()
{
  // NOTE: this requires that Reactive-Transport is done last, or at least
  // after the density of water can be evaluated.  This could be problematic
  // for, e.g., salinity intrusion problems where water density is a function
  // of concentration itself, but should work for all other problems?
  Teuchos::RCP<Epetra_MultiVector> tcc_surf =
    S_->GetW<CompositeVector>(tcc_surf_key_, tag_next_, "state").ViewComponent("cell", true);
  S_->GetEvaluator(mol_dens_surf_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> mol_dens_surf =
    S_->Get<CompositeVector>(mol_dens_surf_key_, tag_next_).ViewComponent("cell", true);

  Teuchos::RCP<Epetra_MultiVector> tcc =
    S_->GetW<CompositeVector>(tcc_key_, tag_next_, "state").ViewComponent("cell", true);
  S_->GetEvaluator(mol_dens_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> mol_dens =
    S_->Get<CompositeVector>(mol_dens_key_, tag_next_).ViewComponent("cell", true);

  int num_aqueous = chemistry_pk_surf_->num_aqueous_components();
  convertConcentrationToAmanzi(*mol_dens_surf, num_aqueous, *tcc_surf, *tcc_surf);
  convertConcentrationToAmanzi(*mol_dens, num_aqueous, *tcc, *tcc);

  chemistry_pk_surf_->set_aqueous_components(tcc_surf);
  chemistry_pk_surf_->Initialize();
  //*tcc_surf = *chemistry_pk_surf_->aqueous_components();

  chemistry_pk_->set_aqueous_components(tcc);
  chemistry_pk_->Initialize();
  //*tcc = *chemistry_pk_->aqueous_components();

  convertConcentrationToATS(*mol_dens_surf, num_aqueous, *tcc_surf, *tcc_surf);
  convertConcentrationToATS(*mol_dens, num_aqueous, *tcc, *tcc);

  transport_pk_surf_->Initialize();
  transport_pk_->Initialize();
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double
MPCCoupledReactiveTransport::get_dt()
{
  double dTtran = coupled_transport_pk_->get_dt();
  double dTchem = coupled_chemistry_pk_->get_dt();

  if (!chem_step_succeeded_ && (dTchem / dTtran > 0.99)) { dTchem *= 0.5; }

  if (dTtran > dTchem) dTtran = dTchem;

  return dTtran;
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool
MPCCoupledReactiveTransport::AdvanceStep(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  chem_step_succeeded_ = false;

  // First we do a transport step.
  bool fail = coupled_transport_pk_->AdvanceStep(t_old, t_new, reinit);
  if (fail) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << coupled_transport_pk_->name() << " failed." << std::endl;
    return fail;
  }

  // Chemistry on the surface
  Teuchos::RCP<Epetra_MultiVector> tcc_surf =
    S_->GetW<CompositeVector>(tcc_surf_key_, tag_next_, "state").ViewComponent("cell", true);
  S_->GetEvaluator(mol_dens_surf_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> mol_dens_surf =
    S_->Get<CompositeVector>(mol_dens_surf_key_, tag_next_).ViewComponent("cell", true);
  fail = advanceChemistry(
    chemistry_pk_surf_, t_old, t_new, reinit, *mol_dens_surf, tcc_surf, *alquimia_surf_timer_);
  changedEvaluatorPrimary(tcc_surf_key_, tag_next_, *S_);
  if (fail) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << chemistry_pk_surf_->name() << " failed." << std::endl;
    return fail;
  } else {
    transport_pk_surf_->debugger()->WriteCellVector("tcc (chem)", *tcc_surf);
    transport_pk_surf_->VV_PrintSoluteExtrema(*tcc_surf, t_new - t_old);
  }

  // Chemistry in the subsurface
  Teuchos::RCP<Epetra_MultiVector> tcc =
    S_->GetW<CompositeVector>(tcc_key_, tag_next_, "state").ViewComponent("cell", true);
  S_->GetEvaluator(mol_dens_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> mol_dens =
    S_->Get<CompositeVector>(mol_dens_key_, tag_next_).ViewComponent("cell", true);
  try {
    fail = advanceChemistry(chemistry_pk_, t_old, t_new, reinit, *mol_dens, tcc, *alquimia_timer_);
    changedEvaluatorPrimary(tcc_key_, tag_next_, *S_);
  } catch (const Errors::Message& chem_error) {
    fail = true;
  }
  if (fail) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << chemistry_pk_->name() << " failed." << std::endl;
    return fail;
  } else {
    transport_pk_->debugger()->WriteCellVector("tcc (chem)", *tcc);
    transport_pk_->VV_PrintSoluteExtrema(*tcc, t_new - t_old);
  }

  chem_step_succeeded_ = true;
  return fail;
};


} // namespace Amanzi
