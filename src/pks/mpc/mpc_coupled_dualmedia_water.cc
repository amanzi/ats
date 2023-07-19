/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "Teuchos_XMLParameterListHelpers.hpp"
#include "EpetraExt_RowMatrixOut.h"

#include "mpc_coupled_dualmedia_water.hh"

namespace Amanzi {

MPCCoupledDualMediaWater::MPCCoupledDualMediaWater(
  Teuchos::ParameterList& FElist,
  const Teuchos::RCP<Teuchos::ParameterList>& plist,
  const Teuchos::RCP<State>& S,
  const Teuchos::RCP<TreeVector>& soln)
  : PK(FElist, plist, S, soln), StrongMPC<PK_BDF_Default>(FElist, plist, S, soln)
{}


void
MPCCoupledDualMediaWater::Setup(const Teuchos::Ptr<State>& S)
{
  // tweak the sub-PK parameter lists
  StrongMPC<PK_BDF_Default>::Setup(S);

  std::cout << "Setup\n";
  Teuchos::Array<std::string> names = plist_->get<Teuchos::Array<std::string>>("PKs order");
  int npks = names.size();
  for (int i = 0; i != npks; ++i) {
    Teuchos::RCP<const CompositeVectorSpace> tmp = solution_->Map().SubVector(i)->Data();
    std::cout << names[i] << " " << tmp << "\n";
  }

  pks_list_->sublist(names[1]).set("coupled to subsurface via head", true);
  // cast the PKs
  integrated_flow_pk_ = Teuchos::rcp_dynamic_cast<StrongMPC<PK_PhysicalBDF_Default>>(sub_pks_[0]);
  AMANZI_ASSERT(integrated_flow_pk_ != Teuchos::null);

  macro_flow_pk_ = sub_pks_[1];
  matrix_flow_pk_ = integrated_flow_pk_->get_subpk(0);
  surf_flow_pk_ = integrated_flow_pk_->get_subpk(1);

  macro_flux_key_ = "macropore-mass_flux";
  matrix_flux_key_ = "water_flux";
}

void
MPCCoupledDualMediaWater::Initialize(const Teuchos::Ptr<State>& S)
{
  // Initialize my timestepper.
  PK_BDF_Default::Initialize(S);
  StrongMPC<PK_BDF_Default>::Initialize(S);

  Comm_ptr_type comm = solution_->Comm();
  auto tvs = Teuchos::rcp(new TreeVectorSpace(comm));

  auto op0 = integrated_flow_pk_->my_operator(Operators::OPERATOR_MATRIX)->Clone();
  auto tvs0 = Teuchos::rcp(new TreeVectorSpace(op0->get_domain_map()));
  tvs->PushBack(tvs0);

  auto op1 = sub_pks_[1]->my_operator(Operators::OPERATOR_MATRIX)->Clone();
  auto tvs1 = Teuchos::rcp(new TreeVectorSpace(op1->get_domain_map()));
  tvs->PushBack(tvs1);

  op_tree_matrix_ = Teuchos::rcp(new Operators::TreeOperator(tvs));

  // we assume that 0 and 1 correspond to matrix and fracture, respectively
  // to avoid modifying original operators, we clone them.

  op_tree_matrix_->set_operator_block(0, 0, op0);
  op_tree_matrix_->set_operator_block(1, 1, op1);

  // off-diagonal blocks are coupled PDEs
  // -- minimum composite vector spaces containing the coupling term
  auto mesh_matrix = S_->GetMesh("domain");
  auto mesh_macropore = S_->GetMesh("macropore");

  auto& mmap0 = solution_->SubVector(0)->SubVector(0)->Data()->ViewComponent("face", false)->Map();
  auto& gmap0 = solution_->SubVector(0)->SubVector(0)->Data()->ViewComponent("face", true)->Map();
  auto& mmap1 = solution_->SubVector(1)->Data()->ViewComponent("face", false)->Map();
  auto& gmap1 = solution_->SubVector(1)->Data()->ViewComponent("face", true)->Map();

  auto cvs_matrix_faces = Teuchos::rcp(new CompositeVectorSpace());
  auto cvs_macropore_faces = Teuchos::rcp(new CompositeVectorSpace());

  cvs_matrix_faces->SetMesh(mesh_matrix)
    ->SetGhosted(true)
    ->AddComponent(
      "face", AmanziMesh::Entity_kind::FACE, Teuchos::rcpFromRef(mmap0), Teuchos::rcpFromRef(gmap0), 1);

  cvs_macropore_faces->SetMesh(mesh_macropore)
    ->SetGhosted(true)
    ->AddComponent(
      "face", AmanziMesh::Entity_kind::FACE, Teuchos::rcpFromRef(mmap1), Teuchos::rcpFromRef(gmap1), 1);

  // create a global problem
  // sub_pks_[0]->my_pde(Operators::PDE_DIFFUSION)->ApplyBCs(true, true, true);
  // sub_pks_[1]->my_pde(Operators::PDE_DIFFUSION)->ApplyBCs(true, true, true);


  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    Teuchos::OSTab tab = vo_->getOSTab();
    *vo_->os() << "matrix:" << std::endl
               << op_tree_matrix_->PrintDiagnostics() << std::endl
               << vo_->color("green") << "Initialization of PK is complete: my dT=" << get_dt()
               << vo_->reset() << std::endl
               << std::endl;
  }
}

void
MPCCoupledDualMediaWater::set_states(const Teuchos::RCP<State>& S,
                                     const Teuchos::RCP<State>& S_inter,
                                     const Teuchos::RCP<State>& S_next)
{
  StrongMPC<PK_BDF_Default>::set_states(S, S_inter, S_next);
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
  Solution_to_State(*u_new, S_next_);

  macro_flow_pk_->FunctionalResidual(
    t_old, t_new, u_old->SubVector(1), u_new->SubVector(1), g->SubVector(1));

  // Evaluate the surface flow residual
  surf_flow_pk_->FunctionalResidual(t_old,
                                    t_new,
                                    u_old->SubVector(0)->SubVector(1),
                                    u_new->SubVector(0)->SubVector(1),
                                    g->SubVector(0)->SubVector(1));


  // All surface to subsurface fluxes have been taken by the subsurface.
  g->SubVector(0)->SubVector(1)->Data()->ViewComponent("cell", false)->PutScalar(0.);
}

// -- Apply preconditioner to u and returns the result in Pu.
int
MPCCoupledDualMediaWater::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                              Teuchos::RCP<TreeVector> Pu)
{
  int ierr;

  return (ierr > 0) ? 0 : 1;
}

// -- Update the preconditioner.
void
MPCCoupledDualMediaWater::UpdatePreconditioner(double t,
                                               Teuchos::RCP<const TreeVector> up,
                                               double h)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon update at t = " << t << std::endl;
}


} // namespace Amanzi
