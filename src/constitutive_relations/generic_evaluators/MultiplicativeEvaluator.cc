/*
  MultiplicativeEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)

*/

#include "MultiplicativeEvaluator.hh"

namespace Amanzi {
namespace Relations {

MultiplicativeEvaluator::MultiplicativeEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  if (!plist.isParameter("evaluator dependencies")) {
    if (plist.isParameter("evaluator dependency suffixes")) {
      const auto& names = plist_.get<Teuchos::Array<std::string> >("evaluator dependency suffixes");
      Key domain = Keys::getDomain(my_key_);
      for (const auto& name : names) {
        dependencies_.insert(Keys::getKey(domain, name));
      }
    } else {
      Errors::Message msg;
      msg << "MultiplicativeEvaluator for: \"" << my_key_ << "\" has no dependencies.";
      Exceptions::amanzi_throw(msg);
    }
  }

  coef_ = plist_.get<double>("coefficient", 1.0);
  positive_ = plist_.get<bool>("enforce positivity", false);
}

Teuchos::RCP<FieldEvaluator>
MultiplicativeEvaluator::Clone() const
{
  return Teuchos::rcp(new MultiplicativeEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void
MultiplicativeEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(dependencies_.size() > 1);
  KeySet::const_iterator key = dependencies_.begin();
  *result = *S->GetPtr<CompositeVector>(*key);
  result->Scale(coef_);
  key++;

  for (auto lcv_name : *result) {
    auto& res_c = *result->ViewComponent(lcv_name, false);
    for (; key!=dependencies_.end(); ++key) {
      const auto& dep_c = *S->Get<CompositeVector>(*key).ViewComponent(lcv_name, false);
      for (int c=0; c!=res_c.MyLength(); ++c) res_c[0][c] *= dep_c[0][c];
    }
    if (positive_) {
      for (int c=0; c!=res_c.MyLength(); ++c) {
        res_c[0][c] = std::max(res_c[0][c], 0.);
      }
    }
  }
}

void
MultiplicativeEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(dependencies_.size() > 1);

  KeySet::const_iterator key = dependencies_.begin();
  while (*key == wrt_key) key++;
  *result = *S->GetPtr<CompositeVector>(*key);
  result->Scale(coef_);
  key++;

  for (; key!=dependencies_.end(); ++key) {
    if (*key != wrt_key) {
      Teuchos::RCP<const CompositeVector> dep = S->GetPtr<CompositeVector>(*key);
      for (CompositeVector::name_iterator lcv=result->begin(); lcv!=result->end(); ++lcv) {
        Epetra_MultiVector& res_c = *result->ViewComponent(*lcv, false);
        const Epetra_MultiVector& dep_c = *dep->ViewComponent(*lcv, false);

        for (int c=0; c!=res_c.MyLength(); ++c) res_c[0][c] *= dep_c[0][c];
      }
    }
  }
}


} // namespace
} // namespace

