/*
  This is the mpc_pk component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling of Transport PK and Chemistry PK.
*/

#include "mpc_reactivetransport.hh"
#include "pk_helpers.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCReactiveTransport::MPCReactiveTransport(Teuchos::ParameterList& pk_tree,
                                          const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                          const Teuchos::RCP<State>& S,
                                          const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_list, S, soln),
    WeakMPC(pk_tree, global_list, S, soln)
{
  chem_step_succeeded_ = true;

  transport_pk_index_ = plist_->get<int>("transport index", 0);
  chemistry_pk_index_ = plist_->get<int>("chemistry index", 1 - transport_pk_index_);

  alquimia_timer_ = Teuchos::TimeMonitor::getNewCounter("alquimia");

  domain_ = plist_->get<std::string>("domain name", "domain");
  tcc_key_ = Keys::readKey(*plist_, domain_, "total component concentration", "total_component_concentration");
  mol_den_key_ = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");
}


void MPCReactiveTransport::Setup()
{
  WeakMPC::Setup();
  cast_sub_pks_();

  S_->Require<CompositeVector,CompositeVectorSpace>(tcc_key_, tag_next_, "state")
    .SetMesh(S_->GetMesh(domain_))->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(tcc_key_, tag_next_);

  S_->Require<CompositeVector,CompositeVectorSpace>(mol_den_key_, tag_next_)
    .SetMesh(S_->GetMesh(domain_))->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireEvaluator(mol_den_key_, tag_next_);
}

void MPCReactiveTransport::cast_sub_pks_()
{
  tranport_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_ATS>(sub_pks_[transport_pk_index_]);
  AMANZI_ASSERT(tranport_pk_ != Teuchos::null);

  chemistry_pk_ = Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(sub_pks_[chemistry_pk_index_]);
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);

  // communicate chemistry engine to transport.
#ifdef ALQUIMIA_ENABLED
  tranport_pk_->SetupAlquimia(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_),
                              chemistry_pk_->chem_engine());
#endif
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MPCReactiveTransport::Initialize()
{
  Teuchos::RCP<Epetra_MultiVector> tcc_copy =
    S_->GetW<CompositeVector>(tcc_key_, tag_next_, "state").ViewComponent("cell", true);

  S_->GetEvaluator(mol_den_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> mol_dens =
    S_->Get<CompositeVector>(mol_den_key_, tag_next_).ViewComponent("cell", true);
  ConvertConcentrationToAmanzi(chemistry_pk_, *mol_dens, *tcc_copy, *tcc_copy);

  chemistry_pk_->Initialize();
  ConvertConcentrationToATS(chemistry_pk_, *mol_dens, *tcc_copy, *tcc_copy);
  tranport_pk_->Initialize();
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCReactiveTransport::get_dt()
{
  double dTtran = tranport_pk_->get_dt();
  double dTchem = chemistry_pk_->get_dt();

  if (!chem_step_succeeded_ && (dTchem / dTtran > 0.99)) {
     dTchem *= 0.5;
  }
  return dTchem;
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool MPCReactiveTransport::AdvanceStep(double t_old, double t_new, bool reinit)
{
  bool fail = false;
  chem_step_succeeded_ = false;

  // First we do a transport step.
  bool pk_fail = tranport_pk_->AdvanceStep(t_old, t_new, reinit);
  if (pk_fail) return pk_fail;

  // Second, we do a chemistry step.
  int num_aq_componets = tranport_pk_->get_num_aqueous_component();

  Teuchos::RCP<Epetra_MultiVector> tcc_copy =
    S_->GetW<CompositeVector>(tcc_key_, tag_next_, "state").ViewComponent("cell", true);

  S_->GetEvaluator(mol_den_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const Epetra_MultiVector> mol_dens =
    S_->Get<CompositeVector>(mol_den_key_, tag_next_).ViewComponent("cell", true);

  pk_fail = AdvanceChemistry(chemistry_pk_, *mol_dens, tcc_copy, t_old, t_new, reinit);
  ChangedEvaluatorPrimary(tcc_key_, tag_next_, *S_);
  if (!pk_fail) chem_step_succeeded_ = true;
  return fail;
};


bool
MPCReactiveTransport::AdvanceChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
        const Epetra_MultiVector& mol_dens,
        Teuchos::RCP<Epetra_MultiVector> tcc_copy,
        double t_old, double t_new, bool reinit)
{
  bool pk_fail = false;
  ConvertConcentrationToAmanzi(chem_pk, mol_dens, *tcc_copy, *tcc_copy);
  chem_pk->set_aqueous_components(tcc_copy);
  {
    auto monitor = Teuchos::rcp(new Teuchos::TimeMonitor(*alquimia_timer_));
    pk_fail = chem_pk->AdvanceStep(t_old, t_new, reinit);
  }
  *tcc_copy = *chem_pk->aqueous_components();
  ConvertConcentrationToATS(chem_pk, mol_dens, *tcc_copy, *tcc_copy);
  return pk_fail;
}


void
MPCReactiveTransport::ConvertConcentrationToAmanzi(
  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
  const Epetra_MultiVector& mol_den,
  const Epetra_MultiVector& tcc_ats,
  Epetra_MultiVector& tcc_amanzi)
{
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S_->GetMesh(chem_pk->domain());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  int num_aq_components = chem_pk->num_aqueous_components();

  // convert from mole fraction[-] to mol/L
  for (int k=0; k<num_aq_components; k++) {
    for (int c=0; c<ncells_owned; c++){
      tcc_amanzi[k][c] = tcc_ats[k][c] * (mol_den[0][c] / 1000.);
    }
  }
}

void
MPCReactiveTransport::ConvertConcentrationToATS(
  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
  const Epetra_MultiVector& mol_den,
  const Epetra_MultiVector& tcc_amanzi,
  Epetra_MultiVector& tcc_ats)
{
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S_->GetMesh(chem_pk->domain());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  int num_aq_components = chem_pk->num_aqueous_components();

  // convert from mole fraction[-] to mol/L
  for (int k=0; k<num_aq_components; k++) {
    for (int c=0; c<ncells_owned; c++){
      tcc_ats[k][c] = tcc_amanzi[k][c] / (mol_den[0][c] / 1000.);
    }
  }
}

}  // namespace Amanzi


