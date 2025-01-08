/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  AdditiveEvaluator is the generic evaluator for adding N other fields.

*/

#include "AdditiveEvaluator.hh"

namespace Amanzi {
namespace Relations {

AdditiveEvaluator::AdditiveEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  if (dependencies_.size() == 0) {
    Errors::Message message;
    message << "AdditiveEvaluator: for " << my_keys_[0].first << " was provided no dependencies";
    throw(message);
  }

  for (const auto& dep : dependencies_) {
    Key variable = dep.first;
    Key variable_tag = Keys::getKey(dep.first, dep.second);
    Key varname = Keys::getVarName(variable);
    if (plist.isParameter(variable_tag + " coefficient")) {
      coefs_[variable_tag] = plist.get<double>(variable_tag + " coefficient");
    } else if (plist.isParameter(variable + " coefficient")) {
      coefs_[variable_tag] = plist.get<double>(variable + " coefficient");
    } else if (plist.isParameter(varname + " coefficient")) {
      coefs_[variable_tag] = plist.get<double>(varname + " coefficient");
    } else {
      coefs_[variable_tag] = 1.0;
    }
  }
  shift_ = plist.get<double>("constant shift", 0.);
  positive_ = plist.get<bool>("enforce positivity", false);
}


Teuchos::RCP<Evaluator>
AdditiveEvaluator::Clone() const
{
  return Teuchos::rcp(new AdditiveEvaluator(*this));
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
AdditiveEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(shift_);

  for (const auto& key_tag : dependencies_) {
    const CompositeVector& dep = S.Get<CompositeVector>(key_tag.first, key_tag.second);
    double coef = coefs_[Keys::getKey(key_tag.first, key_tag.second)];
    result[0]->Update(coef, dep, 1.0);
  }

  if (positive_) {
    for (const auto& name : *result[0]) {
      auto& res = *result[0]->ViewComponent(name, false);
      for (int i = 0; i != res.MyLength(); ++i) res[0][i] = std::max(res[0][i], 0.);
    }
  }
}

void
AdditiveEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(coefs_[Keys::getKey(wrt_key, wrt_tag)]);

  if (positive_) {
    const auto& value = S.Get<CompositeVector>(my_keys_.front().first, my_keys_.front().second);
    for (const auto& name : *result[0]) {
      auto& res = *result[0]->ViewComponent(name, false);
      const auto& value_v = *value.ViewComponent(name, false);
      for (int i = 0; i != res.MyLength(); ++i) {
        if (value_v[0][i] == 0.0) { res[0][i] = 0.; }
      }
    }
  }
}


} // namespace Relations
} // namespace Amanzi
