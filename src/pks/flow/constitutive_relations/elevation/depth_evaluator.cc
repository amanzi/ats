/*
  Evaluates depth of various mesh entities.
    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "depth_model.hh"
#include "depth_evaluator.hh"

namespace Amanzi {
namespace Flow {


DepthEvaluator::DepthEvaluator(Teuchos::ParameterList& plist)
    : IndependentVariableFieldEvaluator(plist)
{}

Teuchos::RCP<FieldEvaluator>
DepthEvaluator::Clone() const {
  return Teuchos::rcp(new DepthEvaluator(*this));
}

// Required methods from IndependentVariableFieldEvaluator
void
DepthEvaluator::UpdateField_(const Teuchos::Ptr<State>& S) {
  if (temporally_variable_ || !computed_once_) {
    CompositeVector& result = *S->GetPtrW<CompositeVector>(my_key_, my_key_);
    for (auto& comp : result) {
      if (comp == "cell") {
        // evaluate depths
        Epetra_MultiVector& depth = *result.ViewComponent("cell",false);
        const AmanziMesh::Mesh& mesh = *result.Mesh();
        DepthModel(mesh, depth);
      } else {
        Errors::Message message;
        message << "DepthEvaluator: Depth components on mesh entities named \"" << comp << "\" are not supported.";
        Exceptions::amanzi_throw(message);
      }
    }
    computed_once_ = true;
  }
}

} //namespace
} //namespace

