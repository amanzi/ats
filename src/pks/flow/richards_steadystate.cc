#include "EpetraExt_RowMatrixOut.h"
#include "richards_steadystate.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_RES_FLAG 0


  RichardsSteadyState::RichardsSteadyState(Teuchos::ParameterList& pk_tree,
                                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                                           const Teuchos::RCP<State>& S,
                                           const Teuchos::RCP<TreeVector>& solution) :
    PK(pk_tree, glist, S, solution),
    Richards(pk_tree, glist, S, solution) {}

void RichardsSteadyState::Setup(const Teuchos::Ptr<State>& S) {
  max_iters_ = plist_->sublist("time integrator").get<int>("max iterations", 10);
  Richards::Setup(S);
}

// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void RichardsSteadyState::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    *vo_->os() << "Precon update at t = " << t << std::endl;
  }

  PK_PhysicalBDF_Default::Solution_to_State(*up, S_next_);
  //PKDefaultBase::solution_to_state(*up, S_next_);

  Teuchos::RCP<const CompositeVector> pres = S_next_ -> GetFieldData(key_);
  // update boundary conditions
  bc_pressure_->Compute(S_next_->time());
  bc_flux_->Compute(S_next_->time());
  UpdateBoundaryConditions_(S_next_.ptr());

  preconditioner_->Init();

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());

  // Create the preconditioner
  Teuchos::RCP<const CompositeVector> rel_perm =
      S_next_->GetFieldData(uw_coef_key_);

  S_next_->GetFieldEvaluator(mass_dens_key_)->HasFieldChanged(S_next_.ptr(), name_);
  Teuchos::RCP<const CompositeVector> rho = S_next_->GetFieldData(mass_dens_key_);
  preconditioner_diff_->SetDensity(rho);

  preconditioner_diff_->SetScalarCoefficient(rel_perm, Teuchos::null);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, pres.ptr());

  // Assemble and precompute the Schur complement for inversion.
  preconditioner_diff_->ApplyBCs(true, true, true);
  
  // // TEST
  // if (S_next_->cycle() == 0 && niter_ == 0) {
  //   // Dump the Schur complement
  //   Teuchos::RCP<Epetra_FECrsMatrix> sc = mfd_preconditioner_->Schur();
  //   std::stringstream filename_s;
  //   filename_s << "schur_" << S_next_->cycle() << ".txt";
  //   EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);

  //   // Dump the Afc^T
  //   Teuchos::RCP<const Epetra_CrsMatrix> Afc = mfd_preconditioner_->Afc();
  //   std::stringstream filename_Afc;
  //   filename_Afc << "Afc_" << S_next_->cycle() << ".txt";
  //   EpetraExt::RowMatrixToMatlabFile(filename_Afc.str().c_str(), *Afc);

  //   std::cout << "CYCLE 0, ITER " << niter_ << "!!!!!!!!" << std::endl;

  //   ChangedSolution();
  //   Teuchos::RCP<TreeVector> up_nc = Teuchos::rcp_const_cast<TreeVector>(up);
  //   Teuchos::RCP<TreeVector> up2 = Teuchos::rcp(new TreeVector(*up));
  //   Teuchos::RCP<TreeVector> f1 = Teuchos::rcp(new TreeVector(*up));
  //   Teuchos::RCP<TreeVector> f2 = Teuchos::rcp(new TreeVector(*up));
  //   FunctionalResidual(S_->time(), S_next_->time(), Teuchos::null, up_nc, f1);

  //   *up2 = *up;
  //   int f = 500;
  //   int c = 99;
  //   (*up_nc->Data())("face",f) =
  //       (*up_nc->Data())("face",f) + 10;
  //   (*up_nc->Data())("cell",c) =
  //       (*up_nc->Data())("cell",c) + 20;
  //   ChangedSolution();
  //   FunctionalResidual(S_->time(), S_next_->time(), Teuchos::null, up_nc, f2);

  //   std::cout.precision(16);
  //   std::cout << "DFDP: " << std::endl;
  //   std::cout << "  L0 = " << (*up2->Data())("face",f);
  //   std::cout << "  L1 = " << (*up_nc->Data())("face",f);
  //   std::cout << "  p0 = " << (*up2->Data())("cell",c);
  //   std::cout << "  p1 = " << (*up_nc->Data())("cell",c) << std::endl;
  //   std::cout << "  fL0 = " << (*f1->Data())("face",f);
  //   std::cout << "  fL1 = " << (*f2->Data())("face",f);
  //   std::cout << "  fp0 = " << (*f1->Data())("cell",c);
  //   std::cout << "  fp1 = " << (*f2->Data())("cell",c) << std::endl;
  //   std::cout << "  DfL = " << (*f2->Data())("face",f) - (*f1->Data())("face",f);
  //   std::cout << "  Dfp = " << (*f2->Data())("cell",c) - (*f1->Data())("cell",c) << std::endl;

  //   double df_dp = ((*f2->Data())("face",f)
  //                   -(*f1->Data())("face",f)) / 10;
  //   std::cout << "DFDP = " << df_dp << std::endl;

  //   // invert
  //   f2->Update(-1., *f1, 1.);
  //   //    f2->Print(std::cout);
  //   ApplyPreconditioner(f2, f1);
  //   std::cout << "  dp = " << 10 << std::endl;
  //   std::cout << "  S^-1 * dfL = " << (*f1->Data())("face",f) << std::endl;
  //   std::cout << "  S^-1 * dfp = " << (*f1->Data())("cell",c) << std::endl;
  //   std::cout << std::endl;
  // }
  
  
  // dump the schur complement
  // Teuchos::RCP<Epetra_FECrsMatrix> sc = mfd_preconditioner_->Schur();
  // std::stringstream filename_s;
  // filename_s << "schur_" << S_next_->cycle() << ".txt";
  // EpetraExt::RowMatrixToMatlabFile(filename_s.str().c_str(), *sc);
  // *vo_->os() << "updated precon " << S_next_->cycle() << std::endl;
  
 
};


// RichardsSteadyState is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void RichardsSteadyState::FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                       Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) {
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  
  double h = t_new - t_old;
  AMANZI_ASSERT(std::abs(S_inter_->time() - t_old) < 1.e-4*h);
  AMANZI_ASSERT(std::abs(S_next_->time() - t_new) < 1.e-4*h);

  // pointer-copy temperature into state and update any auxilary data
  Solution_to_State(*u_new, S_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

  if (dynamic_mesh_) matrix_diff_->SetTensorCoefficient(K_);

#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old
               << " t1 = " << t_new << " h = " << h << std::endl;

  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("p_old"); vnames.push_back("p_new");
  std::vector< Teuchos::Ptr<const CompositeVector> > vecs;
  vecs.push_back(S_inter_->GetFieldData(key_).ptr()); vecs.push_back(u.ptr());
  db_->WriteVectors(vnames, vecs, true);
#endif

  // update boundary conditions
  bc_pressure_->Compute(t_new);
  bc_flux_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

  // evaulate water content, because otherwise it is never done.
  S_next_->GetFieldEvaluator(conserved_key_)->HasFieldChanged(S_next_.ptr(), name_);

#if DEBUG_FLAG
  // dump s_old, s_new
  vnames[0] = "sl_old"; vnames[1] = "sl_new";

  vecs[0] = S_inter_->GetFieldData(Keys::getKey(domain_,"saturation_liquid")).ptr();
  vecs[1] = S_next_->GetFieldData(Keys::getKey(domain_,"saturation_liquid")).ptr();

  if (S_next_->HasField(Keys::getKey(domain_,"saturation_ice"))) {
    vnames.push_back("si_old");
    vnames.push_back("si_new");
    vecs.push_back(S_inter_->GetFieldData(Keys::getKey(domain_,"saturation_ice")).ptr());
    vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"saturation_ice")).ptr());
  }

  vnames.push_back("k_rel");
  vecs.push_back(S_next_->GetFieldData(Keys::getKey(domain_,"relative_permeability")).ptr());

  db_->WriteVectors(vnames,vecs,true);

  db_->WriteVector("res (post diffusion)", res.ptr(), true);
#endif

#if DEBUG_RES_FLAG
  if (niter_ < 23) {
    std::stringstream namestream;
    namestream << "flow_residual_" << niter_;
    *S_next_->GetFieldData(namestream.str(),name_) = *res;

    std::stringstream solnstream;
    solnstream << "flow_solution_" << niter_;
    *S_next_->GetFieldData(solnstream.str(),name_) = *u;
  }
#endif
};

} // namespace
} // namespace
