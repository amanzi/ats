/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! SubgridAggregateEvaluator restricts a field to the subgrid version of the same field.

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

 * `"source domain name`" ``[string]`` Domain name of the source mesh.

ONE OF:
* `"field suffix`" ``[string]`` Set the suffix of the variable
OR
* `"field key`" ``[string]`` **DOMAIN-FIELD_SUFFIX**


 */

#include "SubgridAggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

SubgridAggregateEvaluator::SubgridAggregateEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  domain_ = Keys::getDomain(my_key_);
  source_domain_ = plist_.get<std::string>("source domain name");
  if (Keys::isDomainSet(source_domain_)) { // strip the :*
    source_domain_ = Keys::getDomainSetName(source_domain_);
  }
  var_key_ = Keys::getVarName(my_key_);
  nonlocal_dependencies_ = true; // by definition!
}

Teuchos::RCP<FieldEvaluator>
SubgridAggregateEvaluator::Clone() const
{
  return Teuchos::rcp(new SubgridAggregateEvaluator(*this));
}

// Required methods from SecondaryVariableFieldEvaluator
void
SubgridAggregateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto ds = S->GetDomainSet(source_domain_);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell", false);
  for (const auto& subdomain : *ds) {
    ds->DoImport(subdomain,
                 *S->GetFieldData(Keys::getKey(subdomain, var_key_))->ViewComponent("cell", false),
                 result_v);
  }
}

void
SubgridAggregateEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(1.);
}


void
SubgridAggregateEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  if (dependencies_.size() == 0) {
    auto ds = S->GetDomainSet(source_domain_);
    if (ds->get_referencing_parent() == Teuchos::null) {
      Errors::Message msg;
      msg << "SubgridAggregateEvaluator: DomainSet \"" << source_domain_ << "\" does not have a referencing parent but must have one to aggregate.";
      Exceptions::amanzi_throw(msg);
    }
    if (S->GetMesh(domain_) != ds->get_referencing_parent()) {
      Errors::Message msg;
      msg << "SubgridAggregateEvaluator: DomainSet \"" << source_domain_ << "\" has a referencing parent, but it does not match the aggregate vector's domain, \"" << domain_ << "\"";
      Exceptions::amanzi_throw(msg);
    }

    for (const auto& subdomain : *ds) {
      dependencies_.insert(Keys::getKey(subdomain, var_key_));
    }

    Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);
    my_fac->SetMesh(S->GetMesh(domain_))->SetComponent("cell", AmanziMesh::CELL, 1);

    // check plist for vis or checkpointing control
    bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
    S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
    S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

    // Recurse into the tree to propagate info to leaves.
    for (const auto& subdomain : *ds) {
      auto key = Keys::getKey(subdomain, var_key_);
      S->RequireField(key)->SetMesh(S->GetMesh(subdomain))
        ->AddComponent("cell", AmanziMesh::CELL, 1);
      S->RequireFieldEvaluator(key)->EnsureCompatibility(S);
    }
  }
}


} // namespace
} // namespace

