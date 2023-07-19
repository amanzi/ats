/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! SubgridDisaggregateEvaluator restricts a field to the subgrid version of the same field.
#include "SubgridDisaggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

SubgridDisaggregateEvaluator::SubgridDisaggregateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  domain_index_ = domain;
  domain_set_ = Keys::getDomainSetName(domain);
  source_domain_ = plist_.get<std::string>("source domain name");
  Key var_key = Keys::getVarName(my_keys_.front().first);
  source_key_ = Keys::readKey(plist_, source_domain_, "field", var_key);

  auto tag = my_keys_.front().second;
  dependencies_.insert(KeyTag{ source_key_, tag });
}

Teuchos::RCP<Evaluator>
SubgridDisaggregateEvaluator::Clone() const
{
  return Teuchos::rcp(new SubgridDisaggregateEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
SubgridDisaggregateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  auto ds = S.GetDomainSet(domain_set_);
  ds->doExport(domain_index_,
               *S.Get<CompositeVector>(source_key_, tag).ViewComponent("cell", false),
               *result[0]->ViewComponent("cell", false));
}

void
SubgridDisaggregateEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(1.);
}


// Implements custom EC to use dependencies from subgrid
void
SubgridDisaggregateEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                  my_keys_.front().second);
  if (my_fac.HasComponent("cell")) {
    int num_vectors = my_fac.NumVectors("cell");
    CompositeVectorSpace fac;
    fac.SetMesh(S.GetMesh(source_domain_))->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_vectors);
    EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S, fac);
  }
}


} // namespace Relations
} // namespace Amanzi
