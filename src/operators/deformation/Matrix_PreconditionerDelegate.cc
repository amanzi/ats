/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov) (ATS version)
*/

/*
  License:

  Delegate for setting up and solving a preconditioned system.
*/

#include "errors.hh"
#include "dbc.hh"

#include "Matrix_PreconditionerDelegate.hh"

namespace Amanzi {
namespace Operators {

Matrix_PreconditionerDelegate::Matrix_PreconditionerDelegate(Teuchos::ParameterList& plist)
  : plist_(plist)
{
  InitializeFromPlist_();
}


void
Matrix_PreconditionerDelegate::set_matrix(const Teuchos::RCP<Epetra_RowMatrix>& op)
{
  op_ = op;

  // do this once per op set, as it does some setup work that can be saved?
  if (prec_method_ == TRILINOS_ML) {
    ml_prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*op_, solver_plist_, false));
  }
}


void
Matrix_PreconditionerDelegate::InitializePreconditioner()
{
  if (prec_method_ == TRILINOS_ML) {
    if (ml_prec_->IsPreconditionerComputed() ) ml_prec_->DestroyPreconditioner();
    ml_prec_->SetParameterList(solver_plist_);
    ml_prec_->ComputePreconditioner();

  } else if (prec_method_ == TRILINOS_ILU) {
    ilu_prec_ = Teuchos::rcp(new Ifpack_ILU(&*op_));
    ilu_prec_->SetParameters(solver_plist_);
    ilu_prec_->Initialize();
    ilu_prec_->Compute();

  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    Ifpack factory;
    std::string prectype("ILU");
    int ovl = solver_plist_.get<int>("overlap", 0);
    solver_plist_.set<std::string>("schwarz: combine mode", "Add");
    ifp_prec_ = Teuchos::rcp(factory.Create(prectype, &*op_, ovl));
    ifp_prec_->SetParameters(solver_plist_);
    ifp_prec_->Initialize();
    ifp_prec_->Compute();

#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG) {
    hypre_prec_ = Teuchos::rcp(new Ifpack_Hypre(&*op_));
    Teuchos::RCP<FunctionParameter> functs[8];

    int hypre_coarsen_type = solver_plist_.get<int>("coarsen type", 0);
    functs[0] = Teuchos::rcp(
      new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCoarsenType, hypre_coarsen_type));

    int hypre_verbose = solver_plist_.get<int>("verbosity level", 0);
    functs[1] = Teuchos::rcp(
      new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetPrintLevel, hypre_verbose));

    int hypre_nsmooth = solver_plist_.get<int>("number of smoothing iterations", 2);
    functs[2] = Teuchos::rcp(
      new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetNumSweeps, hypre_nsmooth));

    int hypre_ncycles = solver_plist_.get<int>("number of cycles", 2);
    functs[3] = Teuchos::rcp(
      new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetMaxIter, hypre_ncycles));

    int hypre_relax_type = solver_plist_.get<int>("relax type", 6);
    functs[4] = Teuchos::rcp(
      new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetRelaxType, hypre_relax_type));

    double hypre_strong_threshold = solver_plist_.get<double>("strong threshold", 0.5);
    functs[5] = Teuchos::rcp(new FunctionParameter(
      Preconditioner, &HYPRE_BoomerAMGSetStrongThreshold, hypre_strong_threshold));

    double hypre_tol = solver_plist_.get<double>("tolerance", 1e-12);
    functs[6] =
      Teuchos::rcp(new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetTol, hypre_tol));

    int hypre_cycle_type = solver_plist_.get<int>("cycle type", 1);
    functs[7] = Teuchos::rcp(
      new FunctionParameter(Preconditioner, &HYPRE_BoomerAMGSetCycleType, hypre_cycle_type));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", BoomerAMG);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 8);
    hypre_list.set<Teuchos::RCP<FunctionParameter>*>("Functions", functs);

    hypre_prec_->SetParameters(hypre_list);
    hypre_prec_->Initialize();
    hypre_prec_->Compute();

  } else if (prec_method_ == HYPRE_EUCLID) {
    hypre_prec_ = Teuchos::rcp(new Ifpack_Hypre(&*op_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", Euclid);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);

    hypre_prec_->SetParameters(hypre_list);
    hypre_prec_->Initialize();
    hypre_prec_->Compute();

  } else if (prec_method_ == HYPRE_PARASAILS) {
    hypre_prec_ = Teuchos::rcp(new Ifpack_Hypre(&*op_));

    Teuchos::ParameterList hypre_list;
    hypre_list.set("Preconditioner", ParaSails);
    hypre_list.set("SolveOrPrecondition", Preconditioner);
    hypre_list.set("SetPreconditioner", true);
    hypre_list.set("NumFunctions", 0);

    hypre_prec_->SetParameters(hypre_list);
    hypre_prec_->Initialize();
    hypre_prec_->Compute();
#endif
  }
}


int
Matrix_PreconditionerDelegate::ApplyInverse(const Epetra_MultiVector& b,
                                            Epetra_MultiVector& x) const
{
  int ierr = 0;
  if (prec_method_ == TRILINOS_ML) {
    ierr = ml_prec_->ApplyInverse(b, x);
  } else if (prec_method_ == TRILINOS_ILU) {
    ierr = ilu_prec_->ApplyInverse(b, x);
  } else if (prec_method_ == TRILINOS_BLOCK_ILU) {
    ierr = ifp_prec_->ApplyInverse(b, x);
#ifdef HAVE_HYPRE
  } else if (prec_method_ == HYPRE_AMG || prec_method_ == HYPRE_EUCLID ||
             prec_method_ == HYPRE_PARASAILS) {
    ierr = hypre_prec_->ApplyInverse(b, x);
#endif
  } else {
    ierr = -1;
    AMANZI_ASSERT(0); // this error would have been caught already
  }

  return ierr;
}

void
Matrix_PreconditionerDelegate::InitializeFromPlist_()
{
  // which preconditioner to use
  prec_method_ = PREC_METHOD_NULL;

  std::string precmethodstring = plist_.get<std::string>("preconditioner");

  if (precmethodstring == "ML") {
    prec_method_ = TRILINOS_ML;
    solver_plist_ = plist_.sublist("ML Parameters");

  } else if (precmethodstring == "ILU") {
    prec_method_ = TRILINOS_ILU;
    solver_plist_ = plist_.sublist("ILU Parameters");

  } else if (precmethodstring == "Block ILU") {
    prec_method_ = TRILINOS_BLOCK_ILU;
    solver_plist_ = plist_.sublist("Block ILU Parameters");

#ifdef HAVE_HYPRE
  } else if (precmethodstring == "HYPRE AMG") {
    prec_method_ = HYPRE_AMG;
    solver_plist_ = plist_.sublist("HYPRE AMG Parameters");

  } else if (precmethodstring == "HYPRE Euclid") {
    prec_method_ = HYPRE_EUCLID;
    solver_plist_ = plist_.sublist("HYPRE Euclid Parameters");

  } else if (precmethodstring == "HYPRE ParaSails") {
    prec_method_ = HYPRE_EUCLID;
    solver_plist_ = plist_.sublist("HYPRE ParaSails Parameters");
#endif

  } else {
    // throw an error
#ifdef HAVE_HYPRE
    Errors::Message msg(
      "Matrix_MFD: The specified preconditioner " + precmethodstring +
      " is not supported, we only support ML, ILU, HYPRE AMG, HYPRE Euclid, and HYPRE ParaSails");
#else
    Errors::Message msg("Matrix_MFD: The specified preconditioner " + precmethodstring +
                        " is not supported, we only support ML, and ILU");
#endif
    Exceptions::amanzi_throw(msg);
  }
}


} // namespace Operators
} // namespace Amanzi
