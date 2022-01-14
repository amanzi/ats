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
    coefs_[varname] = plist.get<double>(varname+" coefficient", 1.0);
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

  for (std::map<Key, double>::const_iterator it=coefs_.begin();
       it!=coefs_.end(); ++it) {
    const CompositeVector& dep = S.Get<CompositeVector>(it->first);
    result[0]->Update(it->second, dep, 1.0);
  }
}

void
AdditiveEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(coefs_[wrt_key]);
}


} // namespace
} // namespace

