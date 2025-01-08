/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov) (ATS version)
*/

/*
  License:

  Delegate for setting up and solving a preconditioned system.
*/

#ifndef ATS_MATRIX_PRECONDITIONER_DELEGATE_HH_
#define ATS_MATRIX_PRECONDITIONER_DELEGATE_HH_

#include "Epetra_RowMatrix.h"

#include "ml_MultiLevelPreconditioner.h"
#include "Ifpack.h"
#include "Ifpack_ILU.h"
#include "Ifpack_AdditiveSchwarz.h"


namespace Amanzi {
namespace Operators {

class Matrix_PreconditionerDelegate {
 public:
  // Constructor
  Matrix_PreconditionerDelegate(Teuchos::ParameterList& plist);

  // mutator for the operator
  void set_matrix(const Teuchos::RCP<Epetra_RowMatrix>& op);

  // Accessor for the operator.
  Teuchos::RCP<Epetra_RowMatrix> matrix() { return op_; }

  // Initialize the solver
  void InitializePreconditioner();

  // Apply the inverse.
  int ApplyInverse(const Epetra_MultiVector& b, Epetra_MultiVector& x) const;

 protected:
  void InitializeFromPlist_();

 protected:
  // available solver methods
  enum PrecMethod {
    PREC_METHOD_NULL,
    TRILINOS_ML,
    TRILINOS_ILU,
    TRILINOS_BLOCK_ILU,
    HYPRE_AMG,
    HYPRE_EUCLID,
    HYPRE_PARASAILS
  };
  PrecMethod prec_method_;

  // operator
  Teuchos::RCP<Epetra_RowMatrix> op_;

  // parameters
  Teuchos::ParameterList plist_;
  Teuchos::ParameterList solver_plist_;

  // solver structures
#ifdef HAVE_HYPRE
  Teuchos::RCP<Ifpack_Hypre> hypre_prec_;
#endif
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ml_prec_;
  Teuchos::RCP<Ifpack_ILU> ilu_prec_;
  Teuchos::RCP<Ifpack_Preconditioner> ifp_prec_;
};

} // namespace Operators
} // namespace Amanzi


#endif
