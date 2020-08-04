/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy

  Process kernel for coupling of Transport PK and Chemistry PK.
*/

#include "mpc_coupled_reactive_transport_pk.hh"

namespace Amanzi { 

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
Coupled_ReactiveTransport_PK_ATS::Coupled_ReactiveTransport_PK_ATS(
                                          Teuchos::ParameterList& pk_tree,
                                          const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                                          const Teuchos::RCP<State>& S,
                                          const Teuchos::RCP<TreeVector>& soln) :
  ReactiveTransport_PK_ATS(pk_tree, global_list, S, soln)
 { 

//   storage_created = false;
//   chem_step_succeeded = true;
//   std::string pk_name = pk_tree.name();
  
//   boost::iterator_range<std::string::iterator> res = boost::algorithm::find_last(pk_name,"->"); 
//   if (res.end() - pk_name.end() != 0) boost::algorithm::erase_head(pk_name,  res.end() - pk_name.begin());

// // Create miscaleneous lists.
//   Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(global_list, "PKs", true);
//   crt_pk_list_ = Teuchos::sublist(pk_list, pk_name, true);
  
//   transport_pk_index_ = crt_pk_list_->get<int>("transport index", 0);
//   chemistry_pk_index_ = crt_pk_list_->get<int>("chemistry index", 1 - transport_pk_index_);

//   transport_subcycling_ = crt_pk_list_->get<bool>("transport subcycling", false);


//   /*******************************************************************************************/


//   cast_sub_pks_();

 
  // std::cout<<tranport_pk_subsurface_->name()<<"\n";
  // std::cout<<tranport_pk_overland_->name()<<"\n";
  // std::cout<<chemistry_pk_subsurface_->name()<<"\n";
  // std::cout<<chemistry_pk_overland_->name()<<"\n";

  
 }


void Coupled_ReactiveTransport_PK_ATS::cast_sub_pks_(){

  tranport_pk_ = Teuchos::rcp_dynamic_cast<CoupledTransport_PK>(sub_pks_[transport_pk_index_]);
  AMANZI_ASSERT(tranport_pk_ != Teuchos::null);
  
  chemistry_pk_ = Teuchos::rcp_dynamic_cast<WeakMPC>(sub_pks_[chemistry_pk_index_]);
  AMANZI_ASSERT(chemistry_pk_ != Teuchos::null);

  tranport_pk_subsurface_ = 
    Teuchos::rcp_dynamic_cast<Transport::Transport_PK_ATS>(tranport_pk_->get_subpk(0));
  AMANZI_ASSERT(tranport_pk_subsurface_!= Teuchos::null);
  tranport_pk_overland_ = 
    Teuchos::rcp_dynamic_cast<Transport::Transport_PK_ATS>(tranport_pk_->get_subpk(1));
  AMANZI_ASSERT(tranport_pk_overland_!= Teuchos::null);

  chemistry_pk_subsurface_ = 
    Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(chemistry_pk_->get_subpk(0));
  AMANZI_ASSERT(chemistry_pk_subsurface_!= Teuchos::null);
  chemistry_pk_overland_ = 
    Teuchos::rcp_dynamic_cast<AmanziChemistry::Chemistry_PK>(chemistry_pk_->get_subpk(1));
  AMANZI_ASSERT(chemistry_pk_overland_!= Teuchos::null);

  //std::cout<<tranport_pk_overland_->domain_name()<<" "<<chemistry_pk_overland_->domain_name()<<"\n";
  // std::cout<<tranport_pk_subsurface_->domain_name()<<" "
  AMANZI_ASSERT(tranport_pk_overland_->domain_name() == chemistry_pk_overland_->domain_name());
  AMANZI_ASSERT(tranport_pk_subsurface_->domain_name() == chemistry_pk_subsurface_->domain_name());

}


void Coupled_ReactiveTransport_PK_ATS::Setup(const Teuchos::Ptr<State>& S){

  Amanzi::PK_MPCAdditive<PK>::Setup(S);
  
  cast_sub_pks_();

  // std::cout<<tranport_pk_subsurface_->name()<<"\n";
  // std::cout<<tranport_pk_overland_->name()<<"\n";
  // std::cout<<chemistry_pk_subsurface_->name()<<"\n";
  // std::cout<<chemistry_pk_overland_->name()<<"\n";


  
  // communicate chemistry engine to transport.
#ifdef ALQUIMIA_ENABLED
  tranport_pk_subsurface_->SetupAlquimia(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_subsurface_),
                                            chemistry_pk_subsurface_->chem_engine());
  tranport_pk_overland_->SetupAlquimia(Teuchos::rcp_static_cast<AmanziChemistry::Alquimia_PK>(chemistry_pk_overland_),
                                            chemistry_pk_overland_->chem_engine());
#endif


}


  
// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
double Coupled_ReactiveTransport_PK_ATS::get_dt() {

  dTchem_ = chemistry_pk_->get_dt();
  dTtran_ = tranport_pk_->get_dt();


  if (!chem_step_succeeded && (dTchem_/dTtran_ > 0.99)) {
     dTchem_ *= 0.5;
  } 

  if (dTtran_ > dTchem_) dTtran_= dTchem_; 

  if (transport_subcycling_){ 
    return dTchem_;
  } else {
    return dTtran_;
  }

}

void Coupled_ReactiveTransport_PK_ATS::set_states(const Teuchos::RCP<State>& S,
                                                  const Teuchos::RCP<State>& S_inter,
                                                  const Teuchos::RCP<State>& S_next) {
  //  PKDefaultBase::set_states(S, S_inter, S_next);
  //S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;

  //chemistry_pk_->set_states(S, S_inter, S_next);
  tranport_pk_->set_states(S, S_inter, S_next);

}



void Coupled_ReactiveTransport_PK_ATS::Initialize(const Teuchos::Ptr<State>& S){

  //Amanzi::ReactiveTransport_PK_ATS::Initialize(S);
   
  Key subsurface_domain_key = tranport_pk_subsurface_->domain_name();
  Key overland_domain_key = tranport_pk_overland_->domain_name();
  Key tcc_sub_key = Keys::getKey(subsurface_domain_key, "total_component_concentration");
  Key tcc_over_key = Keys::getKey(overland_domain_key, "total_component_concentration");
  Key sub_mol_den_key = Keys::getKey(subsurface_domain_key,  "molar_density_liquid");
  Key over_mol_den_key = Keys::getKey(overland_domain_key,  "molar_density_liquid");

  Teuchos::RCP<Epetra_MultiVector> tcc_sub = 
    S_->GetFieldData(tcc_sub_key,"state")->ViewComponent("cell", true);
  Teuchos::RCP<const Epetra_MultiVector> mol_den_sub = 
    S_->GetFieldData(sub_mol_den_key)->ViewComponent("cell", true);

  Teuchos::RCP<Epetra_MultiVector> tcc_over = 
    S_->GetFieldData(tcc_over_key,"state")->ViewComponent("cell", true);
  Teuchos::RCP<const Epetra_MultiVector> mol_den_over = 
    S_->GetFieldData(over_mol_den_key)->ViewComponent("cell", true);

  ConvertConcentrationToAmanzi(chemistry_pk_subsurface_, *mol_den_sub,  *tcc_sub,  *tcc_sub);
  ConvertConcentrationToAmanzi(chemistry_pk_overland_, *mol_den_over,  *tcc_over,  *tcc_over);
  
  chemistry_pk_overland_->Initialize(S);
  chemistry_pk_subsurface_->Initialize(S);
  
  ConvertConcentrationToATS(chemistry_pk_subsurface_, *mol_den_sub,  *tcc_sub,  *tcc_sub);
  ConvertConcentrationToATS(chemistry_pk_overland_, *mol_den_over,  *tcc_over,  *tcc_over);
 
  tranport_pk_subsurface_->Initialize(S);
  tranport_pk_overland_->Initialize(S);
  


  // //S_->WriteStatistics(vo_);
  // std::cout<<(*tcc_over)<<"\n"; 
  // exit(0);
    
}


bool Coupled_ReactiveTransport_PK_ATS::AdvanceStep(double t_old, double t_new, bool reinit) {

  Teuchos::OSTab tab = vo_->getOSTab();
    
  bool fail = false;
  chem_step_succeeded = false;
  Key subsurface_domain_key = tranport_pk_subsurface_->domain_name();
  Key overland_domain_key = tranport_pk_overland_->domain_name();
  Key tcc_sub_key = Keys::getKey(subsurface_domain_key, "total_component_concentration");
  Key tcc_over_key = Keys::getKey(overland_domain_key, "total_component_concentration");
  Key sub_mol_den_key = Keys::getKey(subsurface_domain_key,  "molar_density_liquid");
  Key over_mol_den_key = Keys::getKey(overland_domain_key,  "molar_density_liquid");

  double t_initial = S_->initial_time();
  double t_final = S_next_->final_time(); 

  *vo_->os()<< "t_initial "<<t_initial<<" t_final "<<t_final<<"\n";
   
  // if (abs(t_old - t_initial) < 1e-12){
  //   double dt_MPC =  t_final - t_initial;
  //   Teuchos::RCP<Epetra_MultiVector> tcc_sub = 
  //     S_->GetFieldCopyData(tcc_sub_key,"subcycling","state")->ViewComponent("cell", true);
  //   Teuchos::RCP<const Epetra_MultiVector> mol_dens_sub =
  //     S_->GetFieldData(sub_mol_den_key)->ViewComponent("cell", true);

  //   AdvanceChemistry(chemistry_pk_subsurface_, *mol_dens_sub,  tcc_sub,  t_initial, t_initial + 0.5*dt_MPC, reinit);    
  // }
    
  // First we do a transport step.

  bool pk_fail = false;
  pk_fail = tranport_pk_->AdvanceStep(t_old, t_new, reinit);

   
  if (pk_fail){
    Errors::Message message("MPC: Coupled Transport PK returned an unexpected error.");
    Exceptions::amanzi_throw(message);
  }

  Teuchos::RCP<Teuchos::TimeMonitor> local_monitor = Teuchos::rcp(new Teuchos::TimeMonitor(*chem_timer_));
   
  try {
    
    Teuchos::RCP<Epetra_MultiVector> tcc_sub = 
      S_->GetFieldCopyData(tcc_sub_key,"subcycling","state")->ViewComponent("cell", true);
   
    Teuchos::RCP<Epetra_MultiVector> tcc_over =
      S_->GetFieldCopyData(tcc_over_key,"subcycling", "state")->ViewComponent("cell", true);        

    Teuchos::RCP<const Epetra_MultiVector> mol_dens_sub =
      S_->GetFieldData(sub_mol_den_key)->ViewComponent("cell", true);

    Teuchos::RCP<const Epetra_MultiVector> mol_dens_over =
      S_->GetFieldData(over_mol_den_key)->ViewComponent("cell", true);


    AdvanceChemistry(chemistry_pk_subsurface_, *mol_dens_sub,  tcc_sub,  t_old, t_new, reinit);
    AdvanceChemistry(chemistry_pk_overland_,   *mol_dens_over, tcc_over, t_old, t_new, reinit);                    

  }
  catch (const Errors::Message& chem_error) {
    fail = true;
  }

  // if (abs(t_new - t_final) < 1e-12){
  //   double dt_MPC =  t_final - t_initial;
  //   Teuchos::RCP<Epetra_MultiVector> tcc_sub = 
  //     S_->GetFieldCopyData(tcc_sub_key,"subcycling","state")->ViewComponent("cell", true);
  //   Teuchos::RCP<const Epetra_MultiVector> mol_dens_sub =
  //     S_->GetFieldData(sub_mol_den_key)->ViewComponent("cell", true);
    
  //   AdvanceChemistry(chemistry_pk_subsurface_, *mol_dens_sub,  tcc_sub,  t_initial + 0.5*dt_MPC, t_final, reinit);
  // }


  
  local_monitor = Teuchos::null;
  
  return fail;
};


// -----------------------------------------------------------------------------
// 
// -----------------------------------------------------------------------------
void Coupled_ReactiveTransport_PK_ATS::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {

  tranport_pk_->CommitStep(t_old, t_new, S);
  chemistry_pk_->CommitStep(t_old, t_new, S);

}

   

}// namespace
