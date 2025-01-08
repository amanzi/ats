/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  ReciprocalEvaluator is the generic evaluator for dividing two vectors.

*/

#include "ReciprocalEvaluator.hh"

namespace Amanzi {
namespace Relations {

ReciprocalEvaluator::ReciprocalEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  if (dependencies_.size() == 0) {
    Errors::Message message;
    message << "ReciprocalEvaluator: for " << my_keys_[0].first << " was provided no dependencies";
    throw(message);
  }

  if (plist.isParameter("reciprocal")) {
    reciprocal_key_ = plist_.get<std::string>("reciprocal");
  } else {
    Errors::Message msg;
    msg << "ReciprocalEvaluator for: \"" << my_keys_[0].first
        << "\" reciprocal is not defined. No reciprocal parameter.";
    Exceptions::amanzi_throw(msg);
  }

  coef_ = plist_.get<double>("coefficient", 1.0);
  positive_ = plist_.get<bool>("enforce positivity", false);
}

Teuchos::RCP<Evaluator>
ReciprocalEvaluator::Clone() const
{
  return Teuchos::rcp(new ReciprocalEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
ReciprocalEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(dependencies_.size() == 2);

  auto key_tag_numer = dependencies_.begin();
  auto key_tag_denom = dependencies_.begin();
  if (key_tag_numer->first == reciprocal_key_)
    key_tag_numer++;
  else
    key_tag_denom++;

  const auto& numer = S.Get<CompositeVector>(key_tag_numer->first, key_tag_numer->second);
  const auto& denom = S.Get<CompositeVector>(key_tag_denom->first, key_tag_denom->second);

  for (const auto& comp : *result[0]) {
    Epetra_Vector& res_c = *(*result[0]->ViewComponent(comp, false))(0);
    const Epetra_Vector& numer_c = *(*numer.ViewComponent(comp, false))(0);
    const Epetra_Vector& denom_c = *(*denom.ViewComponent(comp, false))(0);
    res_c.ReciprocalMultiply(coef_, denom_c, numer_c, 0);
  }

  if (positive_) {
    for (const auto& comp : *result[0]) {
      Epetra_MultiVector& res_c = *result[0]->ViewComponent(comp, false);
      for (int c = 0; c != res_c.MyLength(); ++c) { res_c[0][c] = std::max(res_c[0][c], 0.); }
    }
  }
}

void
ReciprocalEvaluator::EvaluatePartialDerivative_(const State& S,
                                                const Key& wrt_key,
                                                const Tag& wrt_tag,
                                                const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(dependencies_.size() == 2);

  auto key_tag_numer = dependencies_.begin();
  auto key_tag_denom = dependencies_.begin();
  if (key_tag_numer->first == reciprocal_key_)
    key_tag_numer++;
  else
    key_tag_denom++;

  if (wrt_key == key_tag_denom->first) {
    const auto& numer = S.Get<CompositeVector>(key_tag_numer->first, key_tag_numer->second);
    const auto& denom = S.Get<CompositeVector>(key_tag_denom->first, key_tag_denom->second);

    for (const auto& comp : *result[0]) {
      Epetra_MultiVector& res_c = *result[0]->ViewComponent(comp, false);
      const Epetra_MultiVector& numer_c = *numer.ViewComponent(comp, false);
      const Epetra_MultiVector& denom_c = *denom.ViewComponent(comp, false);
      for (int c = 0; c != res_c.MyLength(); ++c)
        res_c[0][c] = -coef_ * numer_c[0][c] / (denom_c[0][c] * denom_c[0][c]);
    }

  } else if (wrt_key == key_tag_numer->first) {
    const auto& denom = S.Get<CompositeVector>(key_tag_denom->first, key_tag_denom->second);

    for (const auto& comp : *result[0]) {
      Epetra_MultiVector& res_c = *result[0]->ViewComponent(comp, false);
      const Epetra_MultiVector& denom_c = *denom.ViewComponent(comp, false);
      for (int c = 0; c != res_c.MyLength(); ++c) res_c[0][c] = -coef_ / denom_c[0][c];
    }
  }
}


} // namespace Relations
} // namespace Amanzi
