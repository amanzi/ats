/*
  This is the mpc_pk component of the Amanzi code.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling of Transport PK and Chemistry PK.
*/

#include "mpc_reactivetransport_pk.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MPCReactiveTransport::MPCReactiveTransport(Teuchos::ParameterList& pk_tree,
                                          const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                          const Teuchos::RCP<State>& S,
                                          const Teuchos::RCP<TreeVector>& soln)
  : WeakMPC(pk_tree, global_list, S, soln)
{
  chem_step_succeeded_ = true;

  transport_pk_index_ = plist_->get<int>("transport index", 0);
  chemistry_pk_index_ = plist_->get<int>("chemistry index", 1 - transport_pk_index_);

  transport_subcycling_ = plist_->get<bool>("transport subcycling", false);

  domain_ = plist_->get<std::string>("domain name", "domain");
  tcc_key_ = Keys::readKey(*plist_, domain_, "total component concentration", "total_component_concentration");
  mol_den_key_ = Keys::readKey(*plist_, domain_, "molar density liquid", "molar_density_liquid");
}


void MPCReactiveTransport::Setup(const Teuchos::Ptr<State>& S)
{
  Amanzi::PK_MPCAdditive<PK>::Setup(S);
  cast_sub_pks_();

  S->Require<CompositeVector,CompositeVectorSpace>(tcc_key_, Tags::NEXT,  "state").SetMesh(S->GetMesh(domain_))->SetGhosted();
  S->Require<CompositeVector,CompositeVectorSpace>(mol_den_key_, Tags::NEXT).SetMesh(S->GetMesh(domain_))->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S->RequireEvaluator(mol_den_key_);
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
void MPCReactiveTransport::Initialize(const Teuchos::Ptr<State>& S)
{
  Teuchos::RCP<Epetra_MultiVector> tcc_copy =
    S_->GetW<CompositeVector>(tcc_key_,"state").ViewComponent("cell", true);

  S->GetEvaluator(mol_den_key_)->HasFieldChanged(S, name_);
  Teuchos::RCP<const Epetra_MultiVector> mol_dens =
    S_->Get<CompositeVector>(mol_den_key_).ViewComponent("cell", true);

  ConvertConcentrationToAmanzi(chemistry_pk_, *mol_dens, *tcc_copy, *tcc_copy);
  chemistry_pk_->Initialize(S);
  ConvertConcentrationToATS(chemistry_pk_, *mol_dens, *tcc_copy, *tcc_copy);
  tranport_pk_->Initialize(S);
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double MPCReactiveTransport::get_dt()
{
  dTtran_ = tranport_pk_->get_dt();
  dTchem_ = chemistry_pk_->get_dt();

  if (!chem_step_succeeded_ && (dTchem_/dTtran_ > 0.99)) {
     dTchem_ *= 0.5;
  }
  if (dTtran_ > dTchem_) dTtran_= dTchem_;
  return dTchem_;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MPCReactiveTransport::set_dt(double dt)
{
  dTtran_ = dt;
  dTchem_ = dt;
  chemistry_pk_->set_dt(dTchem_);
  tranport_pk_->set_dt(dTtran_);
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
  if (pk_fail){
    Errors::Message message("MPC: Transport PK returned an unexpected error.");
    Exceptions::amanzi_throw(message);
  }

  // Second, we do a chemistry step.
  {
    auto chem_monitor = Teuchos::rcp(new Teuchos::TimeMonitor(*chem_timer_));
    try {
      int num_aq_componets = tranport_pk_ -> num_aqueous_component();

      Teuchos::RCP<Epetra_MultiVector> tcc_copy =
        S_->GetFieldCopyData(tcc_key_,"subcycling","state")->ViewComponent("cell", true);

      S_->GetEvaluator(mol_den_key_)->HasFieldChanged(S_.ptr(), name_);
      Teuchos::RCP<const Epetra_MultiVector> mol_dens =
        S_->Get<CompositeVector>(mol_den_key_).ViewComponent("cell", true);

      pk_fail = AdvanceChemistry(chemistry_pk_, *mol_dens, tcc_copy, t_old, t_new, reinit);

      if (!pk_fail) chem_step_succeeded_ = true;
    }
    catch (const Errors::Message& chem_error) {
      fail = true;
    }
  }
  return fail;
};


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void MPCReactiveTransport::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {

  tranport_pk_->CommitStep(t_old, t_new, S);
  chemistry_pk_->CommitStep(t_old, t_new, S);
}

bool MPCReactiveTransport::AdvanceChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
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


void  MPCReactiveTransport::ConvertConcentrationToAmanzi(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                                                             const Epetra_MultiVector& mol_den,
                                                             const Epetra_MultiVector& tcc_ats,
                                                             Epetra_MultiVector& tcc_amanzi)
{
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S_->GetMesh(chem_pk->domain());
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  int num_aq_components = chem_pk->num_aqueous_components();

  // convert from mole fraction[-] to mol/L
  for (int k=0; k<num_aq_components; k++)
    for (int c=0; c<ncells_owned; c++){
      tcc_amanzi[k][c] = tcc_ats[k][c] * (mol_den[0][c]/ 1000.);
    }
}

void  MPCReactiveTransport::ConvertConcentrationToATS(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
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
      tcc_ats[k][c] = tcc_amanzi[k][c] / ( mol_den[0][c]/ 1000.);
    }
  }
}

}  // namespace Amanzi


