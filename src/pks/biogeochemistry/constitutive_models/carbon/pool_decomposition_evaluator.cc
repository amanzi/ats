/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Pool decomposition, creation of CO2 in soil:
  Koven et al 13, eqn 1, ki*Ci

*/

#include "Epetra_SerialDenseVector.h"

#include "pool_decomposition_evaluator.hh"

namespace Amanzi {
namespace BGC {
namespace BGCRelations {

PoolDecompositionEvaluator::PoolDecompositionEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  carbon_key_ = plist_.get<std::string>("SOM key", "soil_organic_matter");
  dependencies_.insert(carbon_key_);
  decay_key_ = plist_.get<std::string>("pool decay rate key", "soil_carbon_decay_rate");
  dependencies_.insert(decay_key_);

  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("soil co2 key", "soil_co2_production_rate");
  }
}


PoolDecompositionEvaluator::PoolDecompositionEvaluator(const PoolDecompositionEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other),
    carbon_key_(other.carbon_key_),
    decay_key_(other.decay_key_)
{}

Teuchos::RCP<Evaluator>
PoolDecompositionEvaluator::Clone() const
{
  return Teuchos::rcp(new PoolDecompositionEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
PoolDecompositionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                                           const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> carbon_cv = S->GetPtr<CompositeVector>(carbon_key_);
  const AmanziMesh::Mesh& mesh = *carbon_cv->Mesh();

  const Epetra_MultiVector& C = *carbon_cv->ViewComponent("cell", false);
  const Epetra_MultiVector& k =
    *S->GetPtr<CompositeVector>(decay_key_)->ViewComponent("cell", false);
  Epetra_MultiVector& res_c = *result->ViewComponent("cell", false);
  res_c.Multiply(1., C, k, 0.);
}


void
PoolDecompositionEvaluator::EvaluateFieldPartialDerivative_(
  const Teuchos::Ptr<State>& S,
  Key wrt_key,
  const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(0);
}


} // namespace BGCRelations
} // namespace BGC
} // namespace Amanzi
