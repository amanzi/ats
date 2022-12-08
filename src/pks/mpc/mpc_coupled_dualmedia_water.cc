#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "pk_helpers.hh"
#include "mpc_coupled_dualmedia_water.hh"
#include "mpc_surface_subsurface_helpers.hh"

namespace Amanzi {

MPCCoupledDualMediaWater::MPCCoupledDualMediaWater(
  Teuchos::ParameterList& pk_tree,
  const Teuchos::RCP<Teuchos::ParameterList>& plist,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, plist, S, soln),
    StrongMPC<PK_BDF_Default>(pk_tree, plist, S, soln)
{}


void
MPCCoupledDualMediaWater::Setup()
{
  // tweak the sub-PK parameter lists

  macro_flux_key_ = "macropore-water_flux";
  matrix_flux_key_ = "water_flux";
  domain_ss_ = "domain";
  domain_surf_ = "surface";
  domain_macro_ = "macropore";
  
  std::cout << "Setup\n";
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string>>("PKs order");
  int npks = names.size();
  for (int i = 0; i != npks; ++i) {
    Teuchos::RCP<const CompositeVectorSpace> tmp = solution_->Map().SubVector(i)->Data();
  }


  Teuchos::Array<std::string> subnames = pks_list_->sublist(names[0]).get< Teuchos::Array<std::string> >("PKs order");
  
  pks_list_->sublist(subnames[0]).set("coupled to surface via flux", true);
  pks_list_->sublist(subnames[1]).set("coupled to subsurface via flux", true);
  pks_list_->sublist(subnames[1]).set("coupled to macropore via head", true);
  //pks_list_->sublist(names[1]).set("coupled to surface via head", true);
  pks_list_->sublist(names[1]).set("surface domain name", domain_surf_);
  
  ss_flux_key_ =
    Keys::readKey(*plist_, domain_surf_, "surface-subsurface flux",
		  domain_surf_ + "_" + domain_ss_ + "_" + "flux");
  //  if (!S_->HasRecordSet(ss_flux_key_)){
  S_->Require<CompositeVector, CompositeVectorSpace>(ss_flux_key_, tag_next_, ss_flux_key_)
      .SetMesh(S_->GetMesh(domain_surf_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    requireEvaluatorPrimary(ss_flux_key_, tag_next_, *S_);
    //}

  ss_macro_flux_key_ =
    Keys::readKey(*plist_, domain_surf_, "surface-macropore flux",
		  domain_surf_ + "_" + domain_macro_ + "_" + "flux");
  //if (!S_->HasRecordSet(ss_macro_flux_key_)){
  S_->Require<CompositeVector, CompositeVectorSpace>(ss_macro_flux_key_, tag_next_, ss_macro_flux_key_)
      .SetMesh(S_->GetMesh(domain_surf_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    requireEvaluatorPrimary(ss_macro_flux_key_, tag_next_, *S_);
    //}

  StrongMPC<PK_BDF_Default>::Setup();
  
    // cast the PKs
  integrated_flow_pk_ = Teuchos::rcp_dynamic_cast<MPCCoupledWater>(sub_pks_[0]);
  AMANZI_ASSERT(integrated_flow_pk_ != Teuchos::null);

  matrix_flow_pk_ = Teuchos::rcp_dynamic_cast<Flow::Richards>(sub_pks_[1]);
  macro_flow_pk_ = Teuchos::rcp_dynamic_cast<Flow::Richards>(integrated_flow_pk_->get_subpk(0));
  surf_flow_pk_ = Teuchos::rcp_dynamic_cast<Flow::OverlandPressureFlow>(integrated_flow_pk_->get_subpk(1));
  
}

void
MPCCoupledDualMediaWater::Initialize()
{
  // Initialize my timestepper.
  PK_BDF_Default::Initialize();
  StrongMPC<PK_BDF_Default>::Initialize();

  // initialize coupling terms
  S_->GetPtrW<CompositeVector>(ss_flux_key_, tag_next_, ss_flux_key_)->PutScalar(0.);
  S_->GetRecordW(ss_flux_key_, tag_next_, ss_flux_key_).set_initialized();
  changedEvaluatorPrimary(ss_flux_key_, tag_next_, *S_);
  
  S_->GetPtrW<CompositeVector>(ss_macro_flux_key_, tag_next_, ss_macro_flux_key_)->PutScalar(0.);
  S_->GetRecordW(ss_macro_flux_key_, tag_next_, ss_macro_flux_key_).set_initialized();
  changedEvaluatorPrimary(ss_macro_flux_key_, tag_next_, *S_);
  

  Comm_ptr_type comm = solution_->Comm();
  auto tvs = Teuchos::rcp(new TreeVectorSpace(comm));

  

  auto op0 = integrated_flow_pk_->my_operator(Operators::OPERATOR_MATRIX)->Clone();
  // auto tvs0 = Teuchos::rcp(new TreeVectorSpace(op0->get_domain_map()));
  // tvs->PushBack(tvs0);

  // auto op1 = sub_pks_[1]->my_operator(Operators::OPERATOR_MATRIX)->Clone();
  // auto tvs1 = Teuchos::rcp(new TreeVectorSpace(op1->get_domain_map()));
  // tvs->PushBack(tvs1);

  // op_tree_matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  // // we assume that 0 and 1 correspond to matrix and fracture, respectively
  // // to avoid modifying original operators, we clone them.

  // op_tree_matrix_->set_operator_block(0, 0, op0);
  // op_tree_matrix_->set_operator_block(1, 1, op1);

  // // off-diagonal blocks are coupled PDEs
  // // -- minimum composite vector spaces containing the coupling term
  // auto mesh_matrix = S_->GetMesh("domain");
  // auto mesh_macropore = S_->GetMesh("macropore");

  // auto& mmap0 = solution_->SubVector(0)->SubVector(0)->Data()->ViewComponent("face", false)->Map();
  // auto& gmap0 = solution_->SubVector(0)->SubVector(0)->Data()->ViewComponent("face", true)->Map();
  // auto& mmap1 = solution_->SubVector(1)->Data()->ViewComponent("face", false)->Map();
  // auto& gmap1 = solution_->SubVector(1)->Data()->ViewComponent("face", true)->Map();

  // auto cvs_matrix_faces = Teuchos::rcp(new CompositeVectorSpace());
  // auto cvs_macropore_faces = Teuchos::rcp(new CompositeVectorSpace());

  // cvs_matrix_faces->SetMesh(mesh_matrix)
  //   ->SetGhosted(true)
  //   ->AddComponent(
  //     "face", AmanziMesh::FACE, Teuchos::rcpFromRef(mmap0), Teuchos::rcpFromRef(gmap0), 1);

  // cvs_macropore_faces->SetMesh(mesh_macropore)
  //   ->SetGhosted(true)
  //   ->AddComponent(
  //     "face", AmanziMesh::FACE, Teuchos::rcpFromRef(mmap1), Teuchos::rcpFromRef(gmap1), 1);

  // // create a global problem
  // // sub_pks_[0]->my_pde(Operators::PDE_DIFFUSION)->ApplyBCs(true, true, true);
  // // sub_pks_[1]->my_pde(Operators::PDE_DIFFUSION)->ApplyBCs(true, true, true);


  // if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
  //   Teuchos::OSTab tab = vo_->getOSTab();
  //   *vo_->os() << "matrix:" << std::endl
  //              << op_tree_matrix_->PrintDiagnostics() << std::endl
  //              << vo_->color("green") << "Initialization of PK is complete: my dT=" << get_dt()
  //              << vo_->reset() << std::endl
  //              << std::endl;
  // }
}


// -- computes the non-linear functional g = g(t,u,udot)
//    By default this just calls each sub pk FunctionalResidual().
void
MPCCoupledDualMediaWater::FunctionalResidual(double t_old,
                                             double t_new,
                                             Teuchos::RCP<TreeVector> u_old,
                                             Teuchos::RCP<TreeVector> u_new,
                                             Teuchos::RCP<TreeVector> g)
{
  // propagate updated info into state
  Solution_to_State(*u_new, tag_next_);

  const Epetra_MultiVector& matrix_ss_flux = *S_->GetPtr<CompositeVector>(ss_flux_key_, tag_next_)->ViewComponent("cell",false);

  Epetra_MultiVector& source = *S_->GetPtrW<CompositeVector>(ss_macro_flux_key_, tag_next_, ss_macro_flux_key_)->ViewComponent("cell",false);

  std::cout<<"Solution\n";
  std::cout<<"macro cell\n"<<*u_new->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"macro face\n"<<*u_new->SubVector(0)->SubVector(0)->Data()->ViewComponent("face",false)<<"\n";
  std::cout<<"subsurface cell\n"<<*u_new->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"subsurface face\n"<<*u_new->SubVector(1)->Data()->ViewComponent("face",false)<<"\n";
  std::cout<<"surface\n"<<*u_new->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"\n";

  
  matrix_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(1),
                            u_new->SubVector(1), g->SubVector(1));

  std::cout<<"Matrix\n"<<*S_->GetPtr<CompositeVector>(matrix_flux_key_, tag_next_)->ViewComponent("face",false)<<"\n";
  
  CopySubsurfaceToSurface(S_->Get<CompositeVector>(matrix_flux_key_, tag_next_),
			  S_->GetW<CompositeVector>(ss_flux_key_, tag_next_, ss_flux_key_));
  
  std::cout<<"Matrix\n"<<matrix_ss_flux<<"\n";
  std::cout<<"Surface\n"<<source<<"\n";
  
    // Evaluate the surface flow residual
  surf_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(0)->SubVector(1),
                       u_new->SubVector(0)->SubVector(1), g->SubVector(0)->SubVector(1));

  // The residual of the surface flow equation provides the water flux from
  // subsurface to surface.
  
  source = *g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false);
  
  std::cout<<"Surface\n"<<source<<"\n";

  source.Update(-1, matrix_ss_flux, 1);
  changedEvaluatorPrimary(ss_macro_flux_key_, tag_next_, *S_);
  
  // Evaluate the macropore residual, which uses this flux as a Neumann BC.
  macro_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(0)->SubVector(0),
          u_new->SubVector(0)->SubVector(0), g->SubVector(0)->SubVector(0));
  

  // All surface to macropore fluxes have been taken by the macropore
  g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false)->PutScalar(0.);

  std::cout<<"Residual\n";
  std::cout<<"macro cell\n"<<*g->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"macro face\n"<<*g->SubVector(0)->SubVector(0)->Data()->ViewComponent("face",false)<<"\n";
  std::cout<<"subsurfcace cell\n"<<*g->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"subsurface face\n"<<*g->SubVector(1)->Data()->ViewComponent("face",false)<<"\n";
  std::cout<<"surface\n"<<*g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"\n";

}

// -- Apply preconditioner to u and returns the result in Pu.
int
MPCCoupledDualMediaWater::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                              Teuchos::RCP<TreeVector> Pu)
{
  int ierr;
  Pu->PutScalar(0.0);
  //int ok = op_tree_pc_->ApplyInverse(*u, *Pu);
  int ok = StrongMPC<PK_BDF_Default>::ApplyPreconditioner(u, Pu);
  return ok;
}

// -- Update the preconditioner.
void
MPCCoupledDualMediaWater::UpdatePreconditioner(double t,
                                               Teuchos::RCP<const TreeVector> up,
                                               double h)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon update at t = " << t << std::endl;

  StrongMPC<PK_BDF_Default>::UpdatePreconditioner(t, up, h);
  //op_tree_pc_->ComputeInverse();

}


} // namespace Amanzi
