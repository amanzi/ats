/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the porosity after a cell volume change due to mesh deformation.

  Authors: Markus Berndt (berndt@lanl.gov)
*/

#include "porosity_evaluator.hh"

namespace Amanzi {
namespace Deform {
namespace DeformRelations {

PorosityEvaluator::PorosityEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  // add dependency to cell volume
  dependencies_.insert("cell_volume");
  dependencies_.insert("deformation");
}


PorosityEvaluator::PorosityEvaluator(const PorosityEvaluator& other) :
    SecondaryVariableFieldEvaluator(other)
{ }

Teuchos::RCP<FieldEvaluator>
PorosityEvaluator::Clone() const {
  return Teuchos::rcp(new PorosityEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void PorosityEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Epetra_MultiVector& phi_c = *S->GetPtrW<CompositeVector>(my_key_,my_key_)
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& deformation_c = *S->GetPtr<CompositeVector>("deformation")
      ->ViewComponent("cell",false);
  const Epetra_MultiVector& cv = *S->GetPtr<CompositeVector>("cell_volume")
      ->ViewComponent("cell",false);

  // deformation actually stores rock_volume_old
  // new phi = 1 - rock_vol_old / CV_new
  int ncells = phi_c.MyLength();
  for (int c=0; c!=ncells; ++c) {
    phi_c[0][c] = 1. - deformation_c[0][c] / cv[0][c];
    if (phi_c[0][c] < 0.) {
      std::cout << "WARNING: Negative porosity due to too much deformation in cell " << c << std::endl;
      std::cout << "  setting poro = 0" << std::endl;
      phi_c[0][c] = 0.;
    }
  }

}



void PorosityEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  AMANZI_ASSERT(0);
  // not implemented, likely not needed.
}

  
void PorosityEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  // special EnsureCompatibility to add in a evaluator for Porosity
  
  // Call the base class's method since we do not need anything special here
  SecondaryVariableFieldEvaluator::EnsureCompatibility(S);
};


} //namespace
} //namespace
} //namespace

