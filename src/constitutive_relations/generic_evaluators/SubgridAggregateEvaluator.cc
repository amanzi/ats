/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! SubgridAggregateEvaluator restricts a field to the subgrid version of the same field.
#include "SubgridAggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

SubgridAggregateEvaluator::SubgridAggregateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  source_domain_ = plist_.get<std::string>("source domain name");
  if (Keys::isDomainSet(source_domain_)) { // strip the :*
    source_domain_ = Keys::getDomainSetName(source_domain_);
  }
  var_key_ = Keys::getVarName(
    Keys::readKey(plist_, source_domain_, "aggregated", Keys::getVarName(my_keys_.front().first)));
  nonlocal_dependencies_ = true; // by definition!
}

Teuchos::RCP<Evaluator>
SubgridAggregateEvaluator::Clone() const
{
  return Teuchos::rcp(new SubgridAggregateEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
SubgridAggregateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  auto ds = S.GetDomainSet(source_domain_);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  auto dep = dependencies_.begin();
  std::vector<const Epetra_MultiVector*> sources;
  for (const auto& subdomain : *ds) {
    sources.push_back(
      S.Get<CompositeVector>(dep->first, dep->second).ViewComponent("cell", false).get());
    ++dep;
  }
  ds->doImport(sources, result_v);
}

void
SubgridAggregateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                      const Key& wrt_key,
                                                      const Tag& wrt_tag,
                                                      const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(1.);
}


void
SubgridAggregateEvaluator::EnsureEvaluators(State& S)
{
  if (dependencies_.size() == 0) {
    auto ds = S.GetDomainSet(source_domain_);
    Tag dep_tag = Keys::readTag(plist_, my_keys_.front().second);
    if (ds->getReferencingParent() == Teuchos::null) {
      Errors::Message msg;
      msg << "SubgridAggregateEvaluator: DomainSet \"" << source_domain_
          << "\" does not have a referencing parent but must have one to aggregate.";
      Exceptions::amanzi_throw(msg);
    }

    for (const auto& subdomain : *ds) {
      dependencies_.insert(KeyTag{ Keys::getKey(subdomain, var_key_), dep_tag });
    }
  }

  EvaluatorSecondaryMonotypeCV::EnsureEvaluators(S);
}


// Make sure that this vector is set on the referencing parent mesh of the
// domain set.
void
SubgridAggregateEvaluator::EnsureCompatibility_Structure_(State& S)
{
  auto ds = S.GetDomainSet(source_domain_);
  auto& dep_fac = S.Require<CompositeVector, CompositeVectorSpace>(dependencies_.front().first,
                                                                   dependencies_.front().second);
  if (dep_fac.HasComponent("cell")) {
    S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                     my_keys_.front().second)
      .SetMesh(ds->getReferencingParent())
      ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, dep_fac.NumVectors("cell"));
  }

  if (S.GetRecordSet(dependencies_.front().first).subfieldnames()) {
    S.GetRecordSetW(my_keys_.front().first)
      .set_subfieldnames(*S.GetRecordSet(dependencies_.front().first).subfieldnames());
  }
}


// Implements custom EC to use dependencies from subgrid
void
SubgridAggregateEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                               my_keys_.front().second);
  if (fac.HasComponent("cell")) {
    int num_vectors = fac.NumVectors("cell");
    EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S,
                                                              {
                                                                "cell",
                                                              },
                                                              {
                                                                AmanziMesh::Entity_kind::CELL,
                                                              },
                                                              {
                                                                num_vectors,
                                                              });
  }
}


} // namespace Relations
} // namespace Amanzi
