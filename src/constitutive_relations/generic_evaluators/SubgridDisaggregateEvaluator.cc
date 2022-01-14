/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! SubgridDisaggregateEvaluator restricts a field to the subgrid version of the same field.

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "SubgridDisaggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

SubgridDisaggregateEvaluator::SubgridDisaggregateEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  domain_ = Keys::getDomainSetName(my_keys_.front().first);
  source_domain_ = plist_.get<std::string>("source domain name");
  if (Keys::isDomainSet(source_domain_)) { // strip the :*
    source_domain_ = Keys::getDomainSetName(source_domain_);
  }
  var_key_ = Keys::getVarName(my_keys_.front().first);
  source_key_ = Keys::readKey(plist_, source_domain_, "field", var_key_);
}

Teuchos::RCP<Evaluator>
SubgridDisaggregateEvaluator::Clone() const
{
  return Teuchos::rcp(new SubgridDisaggregateEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
SubgridDisaggregateEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  auto ds = S.GetDomainSet(domain_);
  ds->DoExport(domain_,
               *S.Get<CompositeVector>(source_key_).ViewComponent("cell", false),
               *result[0]->ViewComponent("cell", false));
}

void
SubgridDisaggregateEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(1.);
}


void
SubgridDisaggregateEvaluator::EnsureCompatibility(State& S)
{
  if (dependencies_.size() == 0) {
    Key key = my_keys_.front().first;
    Tag tag = my_keys_.front().second;
    dependencies_.insert(KeyTag{source_key_, tag});

    auto ds = S.GetDomainSet(domain_);
    if (ds->get_referencing_parent() == Teuchos::null) {
      Errors::Message msg;
      msg << "SubgridDisaggregateEvaluator: DomainSet \"" << domain_ << "\" does not have a referencing parent but must have one to disaggregate.";
      Exceptions::amanzi_throw(msg);
    }
    if (S.GetMesh(source_domain_) != ds->get_referencing_parent()) {
      Errors::Message msg;
      msg << "SubgridDisaggregateEvaluator: DomainSet \"" << domain_ << "\" has a referencing parent, but it does not match the aggregate vector's domain, \"" << source_domain_ << "\"";
      Exceptions::amanzi_throw(msg);
    }

    S.Require<CompositeVector,CompositeVectorSpace>(key,tag,key)
      .SetMesh(S.GetMesh(Keys::getDomain(key)))
        ->SetComponent("cell", AmanziMesh::CELL, 1);

    S.Require<CompositeVector,CompositeVectorSpace>(source_key_, tag)
      .SetMesh(S.GetMesh(source_domain_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S.RequireEvaluator(source_key_).EnsureCompatibility(S);
  }

  // check plist for vis or checkpointing control
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_Flags_(S);
}


} // namespace
} // namespace

