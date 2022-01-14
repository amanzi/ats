/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! SubgridAggregateEvaluator restricts a field to the subgrid version of the same field.

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "SubgridAggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

SubgridAggregateEvaluator::SubgridAggregateEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomain(my_keys_.front().first);
  source_domain_ = plist_.get<std::string>("source domain name");
  if (Keys::isDomainSet(source_domain_)) { // strip the :*
    source_domain_ = Keys::getDomainSetName(source_domain_);
  }
  var_key_ = Keys::getVarName(my_keys_.front().first);
}

Teuchos::RCP<Evaluator>
SubgridAggregateEvaluator::Clone() const
{
  return Teuchos::rcp(new SubgridAggregateEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
SubgridAggregateEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  auto ds = S.GetDomainSet(source_domain_);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  auto dep = dependencies_.begin();
  std::vector<const Epetra_MultiVector*> sources;
  for (const auto& subdomain : *ds) {
    sources.push_back(S.Get<CompositeVector>(dep->first, dep->second)
                      .ViewComponent("cell", false).get());
    ++dep;
  }
  ds->DoImport(sources, result_v);
}

void
SubgridAggregateEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(1.);
}


void
SubgridAggregateEvaluator::EnsureCompatibility(State& S)
{
  if (dependencies_.size() == 0) {
    auto ds = S.GetDomainSet(source_domain_);
    if (ds->get_referencing_parent() == Teuchos::null) {
      Errors::Message msg;
      msg << "SubgridAggregateEvaluator: DomainSet \"" << source_domain_ << "\" does not have a referencing parent but must have one to aggregate.";
      Exceptions::amanzi_throw(msg);
    }
    if (S.GetMesh(domain_) != ds->get_referencing_parent()) {
      Errors::Message msg;
      msg << "SubgridAggregateEvaluator: DomainSet \"" << source_domain_ << "\" has a referencing parent, but it does not match the aggregate vector's domain, \"" << domain_ << "\"";
      Exceptions::amanzi_throw(msg);
    }

    for (const auto& subdomain : *ds) {
      dependencies_.insert(KeyTag{Keys::getKey(subdomain, var_key_), my_keys_.front().second});
    }

    S.Require<CompositeVector,CompositeVectorSpace>(my_keys_.front().first,
            my_keys_.front().second, my_keys_.front().first)
      .SetMesh(S.GetMesh(domain_))->SetComponent("cell", AmanziMesh::CELL, 1);

    // Recurse into the tree to propagate info to leaves.
    Tag tag = my_keys_.front().second;
    for (const auto& subdomain : *ds) {
      auto key = Keys::getKey(subdomain, var_key_);
      S.Require<CompositeVector,CompositeVectorSpace>(key, tag)
        .SetMesh(S.GetMesh(subdomain))
        ->AddComponent("cell", AmanziMesh::CELL, 1);
      S.RequireEvaluator(key).EnsureCompatibility(S);
    }

    // check plist for vis or checkpointing control
    EvaluatorSecondaryMonotypeCV::EnsureCompatibility_Flags_(S);
  }
}


} // namespace
} // namespace

