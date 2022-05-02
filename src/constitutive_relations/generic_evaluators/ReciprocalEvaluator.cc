/*
  ReciprocalEvaluator is the generic evaluator for dividing two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)

*/

#include "ReciprocalEvaluator.hh"

namespace Amanzi {
namespace Relations {

ReciprocalEvaluator::ReciprocalEvaluator(Teuchos::ParameterList& plist) :
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
      msg << "ReciprocalEvaluator for: \"" << my_key_ << "\" has no dependencies.";
      Exceptions::amanzi_throw(msg);
    }
  }

  if (plist.isParameter("reciprocal")) {
    reciprocal_key_ = plist_.get<std::string>("reciprocal");
  } else {
    Errors::Message msg;
    msg << "ReciprocalEvaluator for: \"" <<my_key_ << "\" reciprocal is not defined. No reciprocal parameter.";
    Exceptions::amanzi_throw(msg);
  }

  coef_ = plist_.get<double>("coefficient", 1.0);
  positive_ = plist_.get<bool>("enforce positivity", false);
  positive_ = plist_.get<bool>("enforce positivity", false);
}

Teuchos::RCP<FieldEvaluator>
ReciprocalEvaluator::Clone() const
{
  return Teuchos::rcp(new ReciprocalEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void
ReciprocalEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(dependencies_.size() == 2);
  result->PutScalar(coef_);

  for (auto key : dependencies_) {
    if (key == reciprocal_key_) {
      for (auto lcv_name : *result) {
        auto& res_c = *result->ViewComponent(lcv_name, false);
        const auto& dep_c = *S->GetFieldData(key)->ViewComponent(lcv_name, false);
        for (int c=0; c!=res_c.MyLength(); ++c) {
          if (std::abs(dep_c[0][c]) < 1e-15) {
            Errors::Message msg;
            msg << "ReciprocalEvaluator for: \"" << my_key_ << "\" has zeros in the denominator. Problem field "<<key<<"\n";
            Exceptions::amanzi_throw(msg);
          }
          res_c[0][c] /= dep_c[0][c];
        }
      }
    } else {
      for (auto lcv_name : *result) {
        auto& res_c = *result->ViewComponent(lcv_name, false);
        const auto& dep_c = *S->GetFieldData(key)->ViewComponent(lcv_name, false);
        for (int c=0; c!=res_c.MyLength(); ++c) {
          res_c[0][c] *= dep_c[0][c];
        }
      }
    }
  }

  if (positive_) {
    for (auto lcv_name : *result) {
      auto& res_c = *result->ViewComponent(lcv_name, false);
      for (int c=0; c!=res_c.MyLength(); ++c) {
        res_c[0][c] = std::max(res_c[0][c], 0.);
      }
    }
  }
}

void
ReciprocalEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  AMANZI_ASSERT(dependencies_.size() == 2);
  KeySet::const_iterator key = dependencies_.begin();

  if ((*key == wrt_key) && (*key!=reciprocal_key_)) {
    key++;
    Teuchos::RCP<const CompositeVector> dep = S->GetFieldData(*key);
    for (CompositeVector::name_iterator lcv=result->begin(); lcv!=result->end(); ++lcv) {
      Epetra_MultiVector& res_c = *result->ViewComponent(*lcv, false);
      const Epetra_MultiVector& dep_c = *dep->ViewComponent(*lcv, false);
      for (int c=0; c!=res_c.MyLength(); ++c) res_c[0][c] = coef_/dep_c[0][c];
    }
  } else {
    *result = *S->GetFieldData(*key);
    result->Scale(coef_);
    key++;
    if ((*key == wrt_key)&&(*key==reciprocal_key_)) {
      Teuchos::RCP<const CompositeVector> dep = S->GetFieldData(*key);
      for (CompositeVector::name_iterator lcv=result->begin(); lcv!=result->end(); ++lcv) {
        Epetra_MultiVector& res_c = *result->ViewComponent(*lcv, false);
        const Epetra_MultiVector& dep_c = *dep->ViewComponent(*lcv, false);
        for (int c=0; c!=res_c.MyLength(); ++c) res_c[0][c] *= -1.0/(dep_c[0][c]*dep_c[0][c]);
      }
    } else {
      result->PutScalar(0.);
    }
  }
}


} // namespace
} // namespace

