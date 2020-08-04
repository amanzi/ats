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
ReactiveTransport_PK_ATS::ReactiveTransport_PK_ATS(Teuchos::ParameterList& pk_tree,
                                          const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                          const Teuchos::RCP<State>& S,
                                          const Teuchos::RCP<TreeVector>& soln) :
  Amanzi::PK_MPCAdditive<PK>(pk_tree, global_list, S, soln)
 { 

  storage_created = false;
  chem_step_succeeded = true;

  std::string pk_name = pk_tree.name();
  
  boost::iterator_range<std::string::iterator> res = boost::algorithm::find_last(pk_name,"->"); 
  if (res.end() - pk_name.end() != 0) boost::algorithm::erase_head(pk_name,  res.end() - pk_name.begin());

// Create miscaleneous lists.
  Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(global_list, "PKs", true);
  rt_pk_list_ = Teuchos::sublist(pk_list, pk_name, true);

  transport_pk_index_ = rt_pk_list_->get<int>("transport index", 0);
  chemistry_pk_index_ = rt_pk_list_->get<int>("chemistry index", 1 - transport_pk_index_);
  
  transport_subcycling_ = rt_pk_list_->get<bool>("transport subcycling", false);

  /*******************************************************************************************/
  

  vo_ = Teuchos::rcp(new Amanzi::VerboseObject(pk_name, *rt_pk_list_));

  // std::cout<<tranport_pk_->domain_name()<<"\n";
  // std::cout<<chemistry_pk_->domain_name()<<"\n";

  // ASSERT(tranport_pk_->domain_name() == chemistry_pk_->domain_name());
  // master_ = 1;  // Transport;
  // slave_ = 0;  // Chemistry;

  chem_timer_ = Teuchos::TimeMonitor::getNewCounter("chemistry");
  alquimia_timer_ = Teuchos::TimeMonitor::getNewCounter("alquimia");
  
}

void ReactiveTransport_PK_ATS::Setup(const Teuchos::Ptr<State>& S){

  Amanzi::PK_MPCAdditive<PK>::Setup(S);
  cast_sub_pks_();
  
}

void ReactiveTransport_PK_ATS::cast_sub_pks_(){

  tranport_pk_ = Teuchos::rcp_dynamic_cast<Transport::Transport_PK_ATS>(sub_pks_[transport_pk_index_]);
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
void ReactiveTransport_PK_ATS::Initialize(const Teuchos::Ptr<State>& S) {

  Teuchos::OSTab tab = vo_->getOSTab();
  //Amanzi::PK_MPCAdditive<PK>::Initialize(S);

  Key domain_name = chemistry_pk_->domain_name();
  Key tcc_key = Keys::getKey(domain_name, "total_component_concentration");
  Key mol_den_key = Keys::getKey(domain_name,  "molar_density_liquid");
  
  Teuchos::RCP<Epetra_MultiVector> tcc_copy =
    S_->GetFieldData(tcc_key,"state")->ViewComponent("cell", true);

  Teuchos::RCP<const Epetra_MultiVector> mol_dens =
    S_->GetFieldData(mol_den_key)->ViewComponent("cell", true);
  
  
  ConvertConcentrationToAmanzi(chemistry_pk_, *mol_dens, *tcc_copy, *tcc_copy);
  
  chemistry_pk_->Initialize(S);
  
  ConvertConcentrationToATS(chemistry_pk_, *mol_dens, *tcc_copy, *tcc_copy);

  tranport_pk_->Initialize(S);


  
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double ReactiveTransport_PK_ATS::get_dt() {

  dTtran_ = tranport_pk_->get_dt();
  dTchem_ = chemistry_pk_->get_dt();

  if (!chem_step_succeeded && (dTchem_/dTtran_ > 0.99)) {
     dTchem_ *= 0.5;
  } 

  if (dTtran_ > dTchem_) dTtran_= dTchem_; 
  
  return dTchem_;
}


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransport_PK_ATS::set_dt(double dt) {
  dTtran_ = dt;
  dTchem_ = dt;
  //dTchem_ = chemistry_pk_ -> get_dt();
  //if (dTchem_ > dTtran_) dTchem_ = dTtran_;
  //if (dTtran_ > dTchem_) dTtran_= dTchem_; 
  chemistry_pk_->set_dt(dTchem_);
  tranport_pk_->set_dt(dTtran_);
}

void ReactiveTransport_PK_ATS::set_states(const Teuchos::RCP<State>& S,
                                          const Teuchos::RCP<State>& S_inter,
                                          const Teuchos::RCP<State>& S_next) {
  //  PKDefaultBase::set_states(S, S_inter, S_next);
  //S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;

  //chemistry_pk_->set_states(S, S_inter, S_next);
  tranport_pk_->set_states(S, S_inter, S_next);

}

// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
bool ReactiveTransport_PK_ATS::AdvanceStep(double t_old, double t_new, bool reinit) {
  bool fail = false;
  chem_step_succeeded = false;

  Key domain_name = tranport_pk_->domain_name();
  Key tcc_key = Keys::getKey(domain_name, "total_component_concentration");
  Key mol_den_key = Keys::getKey(domain_name,  "molar_density_liquid");

  
  // First we do a transport step.
  bool pk_fail = tranport_pk_->AdvanceStep(t_old, t_new, reinit);

  if (pk_fail){
    Errors::Message message("MPC: Transport PK returned an unexpected error.");
    Exceptions::amanzi_throw(message);
  } 

  Teuchos::RCP<Teuchos::TimeMonitor> local_flow_monitor = Teuchos::rcp(new Teuchos::TimeMonitor(*chem_timer_));
  
  // Second, we do a chemistry step.
  try {
    int num_aq_componets = tranport_pk_ -> num_aqueous_component();
    
    Teuchos::RCP<Epetra_MultiVector> tcc_copy =
      S_->GetFieldCopyData(tcc_key,"subcycling","state")->ViewComponent("cell", true);

    Teuchos::RCP<const Epetra_MultiVector> mol_dens =
      S_->GetFieldData(mol_den_key)->ViewComponent("cell", true);
    
    pk_fail = AdvanceChemistry(chemistry_pk_, *mol_dens, tcc_copy, t_old, t_new, reinit);
    
    if (!pk_fail) chem_step_succeeded = true;  
    
    //for (int c=0; c<tcc_copy->MyLength(); c++) std::cout<<(*tcc_copy)[13][c]<<"\n";
    
  }
  catch (const Errors::Message& chem_error) {
    fail = true;
  }

   local_flow_monitor = Teuchos::null;
  
  return fail;
};


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void ReactiveTransport_PK_ATS::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {

  tranport_pk_->CommitStep(t_old, t_new, S);
  chemistry_pk_->CommitStep(t_old, t_new, S);
}
  
bool ReactiveTransport_PK_ATS::AdvanceChemistry(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                                                const Epetra_MultiVector& mol_dens,
                                                Teuchos::RCP<Epetra_MultiVector> tcc_copy,
                                                double t_old, double t_new, bool reinit){

    bool pk_fail = false;


    // convert from mole fraction[-] to mol/L
    // for (int c=0; c<ncells_owned; c++)
    //   for (int k=0; k<num_aq_componets; k++)
    //     (*tcc_copy)[k][c] *= (*mol_dens)[0][c] / 1000.;

    ConvertConcentrationToAmanzi(chem_pk, mol_dens, *tcc_copy, *tcc_copy);

    
    chem_pk->set_aqueous_components(tcc_copy);

    Teuchos::RCP<Teuchos::TimeMonitor> local_flow_monitor = Teuchos::rcp(new Teuchos::TimeMonitor(*alquimia_timer_));
    pk_fail = chem_pk->AdvanceStep(t_old, t_new, reinit);
    local_flow_monitor = Teuchos::null;

    
    *tcc_copy = *chem_pk->aqueous_components();
    
    // convert from mol/L fraction to mole fraction[-]
    // for (int c=0; c<ncells_owned; c++)
    //   for (int k=0; k<num_aq_componets; k++)
    //     (*tcc_copy)[k][c] /= (*mol_dens)[0][c] / 1000.;
    ConvertConcentrationToATS(chem_pk, mol_dens, *tcc_copy, *tcc_copy);

    return pk_fail;

}
  

void  ReactiveTransport_PK_ATS::ConvertConcentrationToAmanzi(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                                                             const Epetra_MultiVector& mol_den,
                                                             const Epetra_MultiVector& tcc_ats,
                                                             Epetra_MultiVector& tcc_amanzi){

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S_->GetMesh(chem_pk->domain_name());
    
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  int num_aq_components = chem_pk->num_aqueous_components();
  
  for (int k=0; k<num_aq_components; k++)
    for (int c=0; c<ncells_owned; c++){
      tcc_amanzi[k][c] = tcc_ats[k][c] * (mol_den[0][c]/ 1000.);
    }
  
}

void  ReactiveTransport_PK_ATS::ConvertConcentrationToATS(Teuchos::RCP<AmanziChemistry::Chemistry_PK> chem_pk,
                                                          const Epetra_MultiVector& mol_den,
                                                          const Epetra_MultiVector& tcc_amanzi,
                                                          Epetra_MultiVector& tcc_ats){

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = S_->GetMesh(chem_pk->domain_name());
    
  int ncells_owned = mesh->num_entities(AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);
  int num_aq_components = chem_pk->num_aqueous_components();
  
  for (int k=0; k<num_aq_components; k++)
    for (int c=0; c<ncells_owned; c++){
      tcc_ats[k][c] = tcc_amanzi[k][c] / ( mol_den[0][c]/ 1000.);
    }

}

}  // namespace Amanzi

 
