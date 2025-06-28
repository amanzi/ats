/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Phong V.V. Le (lepv@ornl.gov)
*/
#include <algorithm>
#include "Key.hh"
#include "Factory.hh"
#include "Function.hh"
#include "qc_relation_multi_evaluator.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {


QCRelationMultiEvaluator::QCRelationMultiEvaluator(Teuchos::ParameterList& plist) : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  cv_key_ = Keys::readKey(plist, domain_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  molar_density_key_ = Keys::readKey(plist, domain_, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ molar_density_key_, tag });
  first_src_key_ = Keys::readKey(plist, domain_, "first source", "first_source");
  dependencies_.insert(KeyTag{ first_src_key_, tag });
  second_src_key_ = Keys::readKey(plist, domain_, "second source", "second_source");
  dependencies_.insert(KeyTag{ second_src_key_, tag });
}

// Required methods from SecondaryVariableMultiEvaluator
void
QCRelationMultiEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result) 
{

  Tag tag = my_keys_.front().second;

  const auto& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const auto& molar_den =
    *S.Get<CompositeVector>(molar_density_key_, tag).ViewComponent("cell", false);
  const auto& first_source =
    *S.Get<CompositeVector>(first_src_key_, tag).ViewComponent("cell", false);
  const auto& second_source =
    *S.Get<CompositeVector>(second_src_key_, tag).ViewComponent("cell", false);    
  
  auto& surf_src = *result[0]->ViewComponent("cell"); // not being reference
  const AmanziMesh::Mesh& mesh = *result[0]->Mesh();
  
  // Loop through each cell
  AmanziMesh::Entity_ID ncells = cv.MyLength();
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {  
    // transport source as a function of discharge and another variable
    surf_src[0][c] = second_source[0][c] * first_source[0][c] * cv[0][c];
  }
}

} // namespace Relations
} // namespace Flow
} // namespace Amanzi 