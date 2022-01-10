/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! SubgridDisaggregateEvaluator restricts a field to the subgrid version of the same field.

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

#include "SubgridDisaggregateEvaluator.hh"

namespace Amanzi {
namespace Relations {

SubgridDisaggregateEvaluator::SubgridDisaggregateEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  domain_ = Keys::getDomainSetName(my_key_);
  source_domain_ = plist_.get<std::string>("source domain name");
  if (Keys::isDomainSet(source_domain_)) { // strip the :*
    source_domain_ = Keys::getDomainSetName(source_domain_);
  }
  source_key_ = Keys::readKey(plist_, source_domain_, "field", Keys::getVarName(my_key_));
}

Teuchos::RCP<FieldEvaluator>
SubgridDisaggregateEvaluator::Clone() const
{
  return Teuchos::rcp(new SubgridDisaggregateEvaluator(*this));
}

// Required methods from SecondaryVariableFieldEvaluator
void
SubgridDisaggregateEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  auto ds = S->GetDomainSet(domain_);
  ds->DoExport(Keys::getDomain(my_key_),
               *S->Get<CompositeVector>(source_key_).ViewComponent("cell", false),
               *result->ViewComponent("cell", false));
}

void
SubgridDisaggregateEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(1.);
}


void
SubgridDisaggregateEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  if (dependencies_.size() == 0) {
    dependencies_.insert(source_key_);

    auto ds = S->GetDomainSet(domain_);
    if (ds->get_referencing_parent() == Teuchos::null) {
      Errors::Message msg;
      msg << "SubgridDisaggregateEvaluator: DomainSet \"" << domain_ << "\" does not have a referencing parent but must have one to disaggregate.";
      Exceptions::amanzi_throw(msg);
    }
    if (S->GetMesh(source_domain_) != ds->get_referencing_parent()) {
      Errors::Message msg;
      msg << "SubgridDisaggregateEvaluator: DomainSet \"" << domain_ << "\" has a referencing parent, but it does not match the aggregate vector's domain, \"" << source_domain_ << "\"";
      Exceptions::amanzi_throw(msg);
    }

    auto my_fac = S->Require<CompositeVector,CompositeVectorSpace>(my_key_, Tags::NEXT,  my_key_);
    my_fac->SetMesh(S->GetMesh(Keys::getDomain(my_key_)))
      ->SetComponent("cell", AmanziMesh::CELL, 1);

    // Check plist for vis or checkpointing control.
    bool io_my_key = plist_.get<bool>("visualize", true);
    S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
    S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

    S->Require<CompositeVector,CompositeVectorSpace>(source_key_, Tags::NEXT)
      ->SetMesh(S->GetMesh(source_domain_))
      ->AddComponent("cell", AmanziMesh::CELL, 1);
    S->RequireFieldEvaluator(source_key_)->EnsureCompatibility(S);
  }
}


} // namespace
} // namespace

