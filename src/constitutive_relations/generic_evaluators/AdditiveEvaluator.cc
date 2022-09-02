/*
  AdditiveEvaluator is the generic evaluator for adding N other fields.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "AdditiveEvaluator.hh"

namespace Amanzi {
namespace Relations {

AdditiveEvaluator::AdditiveEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  if (dependencies_.size() == 0) {
    Errors::Message message;
    message << "AdditiveEvaluator: for " << my_keys_[0].first
                << " was provided no dependencies";
    throw(message);
  }

  for (const auto& dep : dependencies_) {
    Key varname = Keys::getKey(dep.first, dep.second);
    if (plist.isParameter(varname+" coefficient")) {
      coefs_[varname] = plist.get<double>(varname+" coefficient");
    } else if (plist.isParameter(dep.first+" coefficient")) {
      coefs_[varname] = plist.get<double>(dep.first+" coefficient");
    } else {
      coefs_[varname] = 1.0;
    }
  }
  shift_ = plist.get<double>("constant shift", 0.);
}


Teuchos::RCP<Evaluator>
AdditiveEvaluator::Clone() const
{
  return Teuchos::rcp(new AdditiveEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
AdditiveEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(shift_);

  for (const auto& key_tag : dependencies_) {
    const CompositeVector& dep = S.Get<CompositeVector>(key_tag.first, key_tag.second);
    double coef = coefs_[Keys::getKey(key_tag.first, key_tag.second)];
    result[0]->Update(coef, dep, 1.0);
  }
}

void
AdditiveEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(coefs_[Keys::getKey(wrt_key, wrt_tag)]);
}


} // namespace
} // namespace

