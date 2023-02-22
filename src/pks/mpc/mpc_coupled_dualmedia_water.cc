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
  
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string>>("PKs order");
  int npks = names.size();
  for (int i = 0; i != npks; ++i) {
    Teuchos::RCP<const CompositeVectorSpace> tmp = solution_->Map().SubVector(i)->Data();
  }

  Teuchos::Array<std::string> subnames = pks_list_->sublist(names[0]).get< Teuchos::Array<std::string> >("PKs order");
  
  pks_list_->sublist(subnames[0]).set("coupled to surface via flux", true);
  pks_list_->sublist(subnames[1]).set("coupled to subsurface via flux", true);
  pks_list_->sublist(names[1]).set("coupled to surface via head", true);
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

  single_flow_pk_ = Teuchos::rcp_dynamic_cast<Flow::Richards>(sub_pks_[1]);
  embedded_flow_pk_ = Teuchos::rcp_dynamic_cast<Flow::Richards>(integrated_flow_pk_->get_subpk(0));
  surf_flow_pk_ = Teuchos::rcp_dynamic_cast<Flow::OverlandPressureFlow>(integrated_flow_pk_->get_subpk(1));


  domain_mesh_ = S_->GetMesh(domain_ss_);
  surf_mesh_ = S_->GetMesh(domain_surf_);
  macro_mesh_ = S_->GetMesh(domain_macro_);
  
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
  tvs_ = Teuchos::rcp(new TreeVectorSpace(comm));

  op0 = integrated_flow_pk_->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW)->Clone();
  auto tvs0 = Teuchos::rcp(new TreeVectorSpace(op0->get_domain_map()));
  tvs_->PushBack(tvs0);

  // op1 = surf_flow_pk_->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW);
  // auto tvs1 = Teuchos::rcp(new TreeVectorSpace(->get_domain_map()));
  // tvs->PushBack(tvs1);
  
  op2 = single_flow_pk_->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW)->Clone();
  auto tvs2 = Teuchos::rcp(new TreeVectorSpace(op2->get_domain_map()));
  tvs_->PushBack(tvs2);

  op_coupling02 = single_flow_pk_->my_operator(Operators::OPERATOR_PRECONDITIONER_RAW)->Clone();

  GenerateOffDiagonalBlocks();
  
  // // create a global problem
  // // sub_pks_[0]->my_pde(Operators::PDE_DIFFUSION)->ApplyBCs(true, true, true);
  // // sub_pks_[1]->my_pde(Operators::PDE_DIFFUSION)->ApplyBCs(true, true, true);

  solution_->get_map()->Print(std::cout);
  
  tvs_->Print(std::cout);

  auto inv_list = Teuchos::sublist(plist_, "inverse");  
  // create global matrix
  // -- tree matrix (for other MPCs)
  op_tree_pc_ = Teuchos::rcp(new Operators::TreeOperator(tvs_));
  // we assume that 0 and 1 correspond to integrated hydro and matrix, respectively
  // to avoid modifying original operators, we clone them.
  op_tree_pc_->set_operator_block(0, 0, op0);
  //op_tree_pc_->set_operator_block(1, 1, op1);
  op_tree_pc_->set_operator_block(1, 1, op2);
  op_tree_pc_->set_operator_block(0, 1, op_coupling02);
  op_tree_pc_->set_operator_block(1, 0, coupling20_local_op_->global_operator());  
 

  op_tree_pc_->set_inverse_parameters(*inv_list);
  op_tree_pc_->InitializeInverse();

  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "matrix:" << std::endl
               << op_tree_pc_->PrintDiagnostics() << std::endl
               << vo_->color("green") << "Initialization of PK is complete: my dT=" << get_dt()
               << vo_->reset() << std::endl
               << std::endl;
  }

  // op2->set_inverse_parameters(*inv_list);
  // op2->InitializeInverse(); 

}

void MPCCoupledDualMediaWater::GenerateOffDiagonalBlocks(){

  // off-diagonal blocks are coupled PDEs
  // -- minimum composite vector spaces containing the coupling term
  // auto& mmap0_macro_face = solution_->SubVector(0)->SubVector(0)->Data()->ViewComponent("face", false)->Map();
  // auto& gmap0_macro_face = solution_->SubVector(0)->SubVector(0)->Data()->ViewComponent("face", true)->Map();
  // auto& mmap1_mtrx_face = solution_->SubVector(1)->Data()->ViewComponent("face", false)->Map();
  // auto& gmap1_mtrx_face = solution_->SubVector(1)->Data()->ViewComponent("face", true)->Map();

  auto cvs_matrix_faces = Teuchos::rcp(new CompositeVectorSpace());
  auto cvs_macropore_faces = Teuchos::rcp(new CompositeVectorSpace());

  cvs_matrix_faces->SetMesh(domain_mesh_)
    ->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 1);
  cvs_macropore_faces->SetMesh(macro_mesh_)
    ->SetGhosted(true)->AddComponent("face", AmanziMesh::FACE, 1);

  int nsurf_cells = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  std::shared_ptr<std::vector<std::vector<int>>> inds_macro = std::make_shared<std::vector<std::vector<int>>>(nsurf_cells);
  std::shared_ptr<std::vector<std::vector<int>>> inds_matrix = std::make_shared<std::vector<std::vector<int>>>(nsurf_cells);
  std::shared_ptr<std::vector<double>> values = std::make_shared<std::vector<double>>(nsurf_cells);
            
  matrix_local_op_ = *(op2->begin());
  

  std::string name = "Coupling Surface2Matrix: CELL_FACE+CELL";
  coupling02_local_op_ = Teuchos::rcp(new Operators::Op_Cell_FaceCell(name, domain_mesh_));

  for (int sc=0; sc<nsurf_cells; ++sc){
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
    AmanziMesh::Entity_ID_List f_cells;
    domain_mesh_->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &f_cells);
    AMANZI_ASSERT(f_cells.size() == 1);
    AmanziMesh::Entity_ID_List c_faces = domain_mesh_->cell_get_faces(f_cells[0]);
    int num_faces = c_faces.size();
    coupling02_local_op_->matrices[f_cells[0]].Reshape(num_faces+1, num_faces+1);
    coupling02_local_op_->matrices[f_cells[0]] = 0.0;

    (*inds_matrix)[sc].resize(1);
    (*inds_macro)[sc].resize(1);
    (*inds_matrix)[sc][0] = f;
    (*inds_macro)[sc][0] = f;
    (*values)[sc] = -1;
  }

  op_coupling02->OpErase(op_coupling02->begin(), op_coupling02->end());
  op_coupling02->OpPushBack(coupling02_local_op_);
  op_coupling02->SymbolicAssembleMatrix();


  Teuchos::ParameterList oplist;
  coupling20_local_op_ = Teuchos::rcp(new Operators::PDE_CouplingFlux(oplist,
								     cvs_matrix_faces,
								     cvs_macropore_faces,
								     inds_matrix,
								     inds_macro,
								     Teuchos::null));
  
  coupling20_local_op_ -> Setup(values, 1.0);
  coupling20_local_op_ -> UpdateMatrices(Teuchos::null, Teuchos::null);
  
  
  //for (auto op : *op2) std::cout<<op<<"\n"<<op->schema_row()<<" "<<op->schema_col()<<"\n";
  //std::cout<<"\n";
  //for (auto op : *op_coupling02) std::cout<<op<<"\n"<<op->schema_row()<<" "<<op->schema_col()<<"\n";
  

  //op_coupling02->AssembleMatrix();
  //exit(0);


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
  //const Epetra_MultiVector& macro_ss_flux = *S_->GetPtr<CompositeVector>(ss_macro_flux_key_, tag_next_)->ViewComponent("cell",false);

  Epetra_MultiVector& source = *S_->GetPtrW<CompositeVector>(ss_macro_flux_key_, tag_next_, ss_macro_flux_key_)->ViewComponent("cell",false);
  //Epetra_MultiVector& source = *S_->GetPtrW<CompositeVector>(ss_flux_key_, tag_next_, ss_flux_key_)->ViewComponent("cell",false);

  //std::cout<<"Solution\n";
  //std::cout<<"macro cell\n"<<*u_new->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  //Epetra_MultiVector& tmp = *u_new->SubVector(0)->SubVector(0)->Data()->ViewComponent("face",false);
  //std::cout<<"Solution macro face "<<tmp[0][250]<<"\n";
  //std::cout<<"subsurface cell\n"<<*u_new->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  //std::cout<<"Solution subsurface face "<<(*u_new->SubVector(1)->Data()->ViewComponent("face",false))[0][250]<<"\n";
  //std::cout<<"Solution surface "<<(*u_new->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false))[0][0]<<"\n";
  //std::cout<<"\n";
  
  single_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(1),
                            u_new->SubVector(1), g->SubVector(1));

  std::cout<<"single pk cell\n"<<*g->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"single pk face\n" <<*g->SubVector(1)->Data()->ViewComponent("face",false)<<"\n";

  //std::cout<<"Matrix\n"<<*S_->GetPtr<CompositeVector>(matrix_flux_key_, tag_next_)->ViewComponent("face",false)<<"\n";

  //CopySubsurfaceToSurface(S_->Get<CompositeVector>(macro_flux_key_, tag_next_),
  //			  S_->GetW<CompositeVector>(ss_macro_flux_key_, tag_next_, ss_macro_flux_key_));
  
  CopySubsurfaceToSurface(S_->Get<CompositeVector>(matrix_flux_key_, tag_next_),
         		  S_->GetW<CompositeVector>(ss_flux_key_, tag_next_, ss_flux_key_));
  
  //std::cout<<"Matrix\n"<<matrix_ss_flux<<"\n";
  //std::cout<<"Surface\n"<<source<<"\n";
  
    // Evaluate the surface flow residual
  surf_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(0)->SubVector(1),
                       u_new->SubVector(0)->SubVector(1), g->SubVector(0)->SubVector(1));

  // The residual of the surface flow equation provides the water flux from
  // subsurface to surface.
  
  source = *g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false);
  
  //std::cout<<"Surface\n"<<source<<"\n";

  source.Update(-1, matrix_ss_flux, 1);
  //source.Update(-1, macro_ss_flux, 1);
  
  changedEvaluatorPrimary(ss_macro_flux_key_, tag_next_, *S_);
  //changedEvaluatorPrimary(ss_flux_key_, tag_next_, *S_);
  
  // Evaluate the macropore residual, which uses this flux as a Neumann BC.
  embedded_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(0)->SubVector(0),
          u_new->SubVector(0)->SubVector(0), g->SubVector(0)->SubVector(0));
  

  // All surface to macropore fluxes have been taken by the macropore
  g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false)->PutScalar(0.);

  //std::cout<<"Residual\n";
  //std::cout<<"macro cell\n"<<*g->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"Residual macro face "<<(*g->SubVector(0)->SubVector(0)->Data()->ViewComponent("face",false))[0][250]<<"\n";
  //std::cout<<"subsurfcace cell\n"<<*g->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"Residual subsurface face " <<(*g->SubVector(1)->Data()->ViewComponent("face",false))[0][250]<<"\n";
  std::cout<<"Residual surface "   <<(*g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false))[0][0]<<"\n";
  std::cout<<"\n";



  //   // propagate updated info into state
  // Solution_to_State(*u_new, tag_next_);

  // u_old->get_map()->Print(std::cout);

  // // const Epetra_MultiVector& matrix_ss_flux = *S_->GetPtr<CompositeVector>(ss_flux_key_, tag_next_)->ViewComponent("cell",false);
  // const Epetra_MultiVector& macro_ss_flux = *S_->GetPtr<CompositeVector>(ss_macro_flux_key_, tag_next_)->ViewComponent("cell",false);

  // // Epetra_MultiVector& source = *S_->GetPtrW<CompositeVector>(ss_macro_flux_key_, tag_next_, ss_macro_flux_key_)->ViewComponent("cell",false);
  // Epetra_MultiVector& source = *S_->GetPtrW<CompositeVector>(ss_flux_key_, tag_next_, ss_flux_key_)->ViewComponent("cell",false);

  // //std::cout<<"Solution\n";
  // //std::cout<<"macro cell\n"<<*u_new->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  // //Epetra_MultiVector& tmp = *u_new->SubVector(0)->SubVector(0)->Data()->ViewComponent("face",false);
  // //std::cout<<"Solution macro face "<<tmp[0][250]<<"\n";
  // //std::cout<<"subsurface cell\n"<<*u_new->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  // //std::cout<<"Solution subsurface face "<<(*u_new->SubVector(1)->Data()->ViewComponent("face",false))[0][250]<<"\n";
  // //std::cout<<"Solution surface "<<(*u_new->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false))[0][0]<<"\n";
  // //std::cout<<"\n";
  
  // single_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(1),
  //                           u_new->SubVector(1), g->SubVector(1));

  // std::cout<<"single pk cell\n"<<*g->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  // std::cout<<"single pk face\n" <<*g->SubVector(1)->Data()->ViewComponent("face",false)<<"\n";

  // //std::cout<<"Matrix\n"<<*S_->GetPtr<CompositeVector>(matrix_flux_key_, tag_next_)->ViewComponent("face",false)<<"\n";

  // CopySubsurfaceToSurface(S_->Get<CompositeVector>(macro_flux_key_, tag_next_),
  // 			  S_->GetW<CompositeVector>(ss_macro_flux_key_, tag_next_, ss_macro_flux_key_));
  
  // // CopySubsurfaceToSurface(S_->Get<CompositeVector>(matrix_flux_key_, tag_next_),
  // //       		  S_->GetW<CompositeVector>(ss_flux_key_, tag_next_, ss_flux_key_));
  
  // //std::cout<<"Matrix\n"<<matrix_ss_flux<<"\n";
  // //std::cout<<"Surface\n"<<source<<"\n";
  
  //   // Evaluate the surface flow residual
  // surf_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(0)->SubVector(1),
  //                      u_new->SubVector(0)->SubVector(1), g->SubVector(0)->SubVector(1));

  // // The residual of the surface flow equation provides the water flux from
  // // subsurface to surface.
  
  // source = *g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false);
  
  // //std::cout<<"Surface\n"<<source<<"\n";

  // // source.Update(-1, matrix_ss_flux, 1);
  // source.Update(-1, macro_ss_flux, 1);
  
  // //changedEvaluatorPrimary(ss_macro_flux_key_, tag_next_, *S_);
  // changedEvaluatorPrimary(ss_flux_key_, tag_next_, *S_);
  
  // // Evaluate the macropore residual, which uses this flux as a Neumann BC.
  // embedded_flow_pk_->FunctionalResidual(t_old, t_new, u_old->SubVector(0)->SubVector(0),
  //         u_new->SubVector(0)->SubVector(0), g->SubVector(0)->SubVector(0));
  

  // // All surface to macropore fluxes have been taken by the macropore
  // g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false)->PutScalar(0.);

  // //std::cout<<"Residual\n";
  // //std::cout<<"macro cell\n"<<*g->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  // std::cout<<"Residual macro face "<<(*g->SubVector(0)->SubVector(0)->Data()->ViewComponent("face",false))[0][250]<<"\n";
  // //std::cout<<"subsurfcace cell\n"<<*g->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  // std::cout<<"Residual subsurface face " <<(*g->SubVector(1)->Data()->ViewComponent("face",false))[0][250]<<"\n";
  // std::cout<<"Residual surface "   <<(*g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false))[0][0]<<"\n";
  // std::cout<<"\n";

}

// -- Apply preconditioner to u and returns the result in Pu.
int
MPCCoupledDualMediaWater::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                              Teuchos::RCP<TreeVector> Pu)
{
  int ierr;

  std::cout<<"Vector u\n";
  std::cout<<"embeded cell\n"<<*u->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"emdbbed face\n"<<*u->SubVector(0)->SubVector(0)->Data()->ViewComponent("face",false)<<"\n";
  std::cout<<"surface\n"   <<*u->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"single cell\n"<<*u->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"singleface face\n" <<*u->SubVector(1)->Data()->ViewComponent("face",false)<<"\n";
  std::cout<<"\n";

  Comm_ptr_type comm = solution_->Comm();
  Teuchos::RCP<TreeVector> u_reduced = Teuchos::rcp(new TreeVector(tvs_));
  Teuchos::RCP<TreeVector>Pu_reduced = Teuchos::rcp(new TreeVector(comm));
    
  *u_reduced->SubVector(0)->Data() = *u->SubVector(0)->SubVector(0)->Data();
  *u_reduced->SubVector(1)->Data() = *u->SubVector(1)->Data();


  // std::cout<<"Vector u_reduced\n";
  // std::cout<<"macro cell\n"<<*u_reduced->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  // std::cout<<"macro face\n"<<*u_reduced->SubVector(0)->Data()->ViewComponent("face",false)<<"\n";
  // std::cout<<"subsurfcace cell\n"<<*u_reduced->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  // std::cout<<"subsurface face\n" <<*u_reduced->SubVector(1)->Data()->ViewComponent("face",false)<<"\n";
  // std::cout<<"\n";


  // u_reduced->PushBack(u->SubVector(1));
  Pu_reduced->PushBack(Pu->SubVector(0)->SubVector(0));
  Pu_reduced->PushBack(Pu->SubVector(1));

  Pu->PutScalar(0.0);

  Pu_reduced->PutScalar(0.0);
  int ok = op_tree_pc_->ApplyInverse(*u_reduced, *Pu_reduced);

  //int ok = StrongMPC<PK_BDF_Default>::ApplyPreconditioner(u, Pu);

  const Epetra_MultiVector& subsurface_vec_f = *Pu->SubVector(1)->Data()->ViewComponent("face",false); 
  Epetra_MultiVector& surface_vec = *Pu->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false);
  unsigned int ncells_surface = surface_vec.MyLength();

  for (unsigned int c = 0; c != ncells_surface; ++c) {
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
#ifdef ENABLE_DBC
    AmanziMesh::Entity_ID_List cells;
    domain_mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    AMANZI_ASSERT(cells.size() == 1);
#endif
    surface_vec[0][c] = subsurface_vec_f[0][f];
  }

  // std::cout<<"Vector Pu_reduced\n";
  // std::cout<<"macro cell\n"<<*Pu_reduced->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  // std::cout<<"macro face\n"<<*Pu_reduced->SubVector(0)->Data()->ViewComponent("face",false)<<"\n";
  // std::cout<<"subsurfcace cell\n"<<*Pu_reduced->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  // std::cout<<"subsurface face\n" <<*Pu_reduced->SubVector(1)->Data()->ViewComponent("face",false)<<"\n";
  // std::cout<<"\n";

  
  std::cout<<"Vector Pu\n";
  std::cout<<"macro cell\n"<<*Pu->SubVector(0)->SubVector(0)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"macro face\n"<<*Pu->SubVector(0)->SubVector(0)->Data()->ViewComponent("face",false)<<"\n";
  std::cout<<"surface\n"   <<*Pu->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"subsurfcace cell\n"<<*Pu->SubVector(1)->Data()->ViewComponent("cell",false)<<"\n";
  std::cout<<"subsurface face\n" <<*Pu->SubVector(1)->Data()->ViewComponent("face",false)<<"\n";
  std::cout<<"\n";


  //  exit(-1);

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
  UpdateOffDiagonalBlocks();
  
  op_tree_pc_->ComputeInverse();

  //exit(0); 
}

void
MPCCoupledDualMediaWater::UpdateOffDiagonalBlocks(){

  std::cout<<"UpdateOffDiagonalBlocks\n";
  int ncells_owned = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  for (int sc=0; sc<ncells_owned; ++sc){    
    AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, sc);
    std::cout<<"FACE: "<<f<<"\n";
    AmanziMesh::Entity_ID_List f_cells;
    domain_mesh_->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &f_cells);
    AMANZI_ASSERT(f_cells.size() == 1);
    std::cout << matrix_local_op_->matrices[f_cells[0]]<<"\n";
    AmanziMesh::Entity_ID_List c_faces = domain_mesh_->cell_get_faces(f_cells[0]);
    int num_faces = c_faces.size();
    for (int i=0; i<num_faces; ++i){
      if (c_faces[i]==f){
	for (int j=0; j<num_faces+1; j++){
	  (coupling02_local_op_->matrices[f_cells[0]])(i,j) = (matrix_local_op_->matrices_shadow[f_cells[0]])(i,j);
	  //(coupling02_local_op_->matrices[f_cells[0]])(i,j) = 88800.0;
	  if (j==i){
	    (matrix_local_op_->matrices[f_cells[0]])(i,j) = 1.;
	  }else{
	    (matrix_local_op_->matrices[f_cells[0]])(i,j) = 0;
	  }
	}
      }
    }
    std::cout << matrix_local_op_->matrices[f_cells[0]] <<"\n";
    std::cout << coupling02_local_op_->matrices[f_cells[0]] <<"\n";
    std::cout << matrix_local_op_->matrices_shadow[f_cells[0]] <<"\n";
    // std::cout << coupling02_local_op_->matrices_shadow[f_cells[0]] <<"\n";
    
  }
 
  //op_coupling02->AssembleMatrix();

  
}



} // namespace Amanzi
