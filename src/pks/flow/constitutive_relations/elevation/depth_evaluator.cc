/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "depth_model.hh"
#include "depth_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {


DepthEvaluator::DepthEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorIndependentCV(plist)
{
  algorithm_ = plist_.get<std::string>("algorithm", "mean face centroid");
  if (!(algorithm_ == "mean face centroid" || algorithm_ == "cell centroid")) {
    Errors::Message msg;
    msg << "In evaluator DepthEvaluator for \"" << my_key_ << "\": invalid algorithm \""
        << algorithm_ << "\", valid are \"mean face centroid\" and \"cell centroid\"";
    Exceptions::amanzi_throw(msg);
  }
}

Teuchos::RCP<Evaluator>
DepthEvaluator::Clone() const
{
  return Teuchos::rcp(new DepthEvaluator(*this));
}

// Required methods from IndependentVariableEvaluator
void
DepthEvaluator::Update_(State& S)
{
  if (temporally_variable_ || !computed_once_) {
    CompositeVector& result = S.GetW<CompositeVector>(my_key_, my_tag_, my_key_);
    for (auto& comp : result) {
      if (comp == "cell") {
        // evaluate depths
        Epetra_MultiVector& depth = *result.ViewComponent("cell", false);
        const AmanziMesh::Mesh& mesh = *result.Mesh();
        if (algorithm_ == "mean face centroid") {
          computeDepth_MeanFaceCentroid(mesh, depth);
        } else {
          computeDepth_CellCentroid(mesh, depth);
        }
      } else {
        Errors::Message message;
        message << "DepthEvaluator: Depth components on mesh entities named \"" << comp
                << "\" are not supported.";
        Exceptions::amanzi_throw(message);
      }
    }
    computed_once_ = true;
  }
}

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
