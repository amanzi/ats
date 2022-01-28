/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MultiplicativeEvaluator.hh"

namespace Amanzi {
namespace Relations {

MultiplicativeEvaluator::MultiplicativeEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  if (dependencies_.size() == 0) {
    Errors::Message message;
    message << "MultiplicativeEvaluator: for " << my_keys_[0].first
                << " was provided no dependencies";
    throw(message);
  }

  coef_ = plist_.get<double>("coefficient", 1.0);
  positive_ = plist_.get<bool>("enforce positivity", false);
}


Teuchos::RCP<Evaluator>
MultiplicativeEvaluator::Clone() const
{
  return Teuchos::rcp(new MultiplicativeEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
MultiplicativeEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(dependencies_.size() > 1);
  auto key_tag = dependencies_.begin();
  *result[0] = S.Get<CompositeVector>(key_tag->first, key_tag->second);
  result[0]->Scale(coef_);
  key_tag++;

  for (const auto& lcv_name : *result[0]) {
    auto& res_c = *result[0]->ViewComponent(lcv_name, false);
    for (; key_tag!=dependencies_.end(); ++key_tag) {
      const Epetra_MultiVector& dep_c = *S.Get<CompositeVector>(key_tag->first, key_tag->second)
        .ViewComponent(lcv_name, false);
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
MultiplicativeEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag,
        const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(dependencies_.size() > 1);

  auto key_tag = dependencies_.begin();
  KeyTag wrt(wrt_key, wrt_tag);
  while (*key_tag == wrt) key_tag++;
  *result[0] = S.Get<CompositeVector>(key_tag->first, key_tag->second);
  result[0]->Scale(coef_);
  key_tag++;

  for (; key_tag!=dependencies_.end(); ++key_tag) {
    if (*key_tag != wrt) {
      const CompositeVector& dep = S.Get<CompositeVector>(key_tag->first, key_tag->second);
      for (const auto& comp : *result[0]) {
        Epetra_MultiVector& res_c = *result[0]->ViewComponent(comp, false);
        const Epetra_MultiVector& dep_c = *dep.ViewComponent(comp, false);
        // for (int c=0; c!=res_c.MyLength(); ++c) res_c[0][c] *= dep_c[0][c];
        res_c.Multiply(1, res_c, dep_c, 0);
      }
    }
  }
}


} // namespace
} // namespace

