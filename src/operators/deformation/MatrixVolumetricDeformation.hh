/*
  License:
  Authors: Ethan Coon (ecoon@lanl.gov) (ATS version)

  Deformation optimization matrix
*/

#ifndef OPERATORS_MATRIX_VOLUMETRIC_DEFORMATION_HH_
#define OPERATORS_MATRIX_VOLUMETRIC_DEFORMATION_HH_


#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_CrsMatrix.h"

#include "Mesh.hh"
#include "CompositeVector.hh"
#include "CompositeMatrix.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace Operators {

class MatrixVolumetricDeformation : public CompositeMatrix {

 public:

  MatrixVolumetricDeformation(Teuchos::ParameterList& plist,
          const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);

  MatrixVolumetricDeformation(const MatrixVolumetricDeformation& other);


  // Vector space of the Matrix's domain.
  virtual const CompositeVectorSpace& DomainMap() const {
    return *domain_; }

  // Vector space of the Matrix's range.
  virtual const CompositeVectorSpace& RangeMap() const {
    return *range_; }

  // Virtual copy constructor.
  virtual Teuchos::RCP<CompositeMatrix> Clone() const {
    return Teuchos::rcp(new MatrixVolumetricDeformation(*this)); }

  // Apply matrix, b <-- Ax
  virtual int Apply(const CompositeVector& x,
                     CompositeVector& b) const;

  // Apply the inverse, x <-- A^-1 b
  virtual int ApplyInverse(const CompositeVector& b,
                            CompositeVector& x) const;

  // This is a Normal equation, so we need to apply N^T to the rhs
  void ApplyRHS(const CompositeVector& x_cell,
                const Teuchos::Ptr<CompositeVector>& x_node,
                const Teuchos::Ptr<const AmanziMesh::Entity_ID_List>& fixed_nodes) const;

  // Sets up the solver.
  void set_inverse_parameters();
  void Assemble(const Teuchos::Ptr<const AmanziMesh::Entity_ID_List>& fixed_nodes);

 protected:
  void InitializeFromOptions_();
  void PreAssemble_();

 protected:
  // local data
  Teuchos::RCP<CompositeVectorSpace> range_;
  Teuchos::RCP<CompositeVectorSpace> domain_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList plist_;
  Teuchos::RCP<Epetra_CrsMatrix> operator_;
  Teuchos::RCP<Epetra_CrsMatrix> operatorPre_;
  Teuchos::RCP<Epetra_CrsMatrix> dVdz_;
  int max_nnode_neighbors_;

  // parameters to control optimization
  double diagonal_shift_;
  double smoothing_;
};


}
}



#endif
