/*
  AdditiveEvaluator is the generic evaluator for adding N other fields.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "AdditiveEvaluator.hh"

namespace Amanzi {
namespace Relations {

AdditiveEvaluator::AdditiveEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondary(plist)
{
  Teuchos::Array<std::string> deps =
      plist_.get<Teuchos::Array<std::string> >("evaluator dependencies");

  for (auto& dep : deps) {
    Key pname = dep + std::string(" coefficient");
    coefs_[dep] = plist.get<double>(pname, 1.0);
  }  
}


AdditiveEvaluator::AdditiveEvaluator(const AdditiveEvaluator& other) :
    EvaluatorSecondary(other),
    coefs_(other.coefs_) {}

Teuchos::RCP<Evaluator>
AdditiveEvaluator::Clone() const
{
  return Teuchos::rcp(new AdditiveEvaluator(*this));
}

// Required methods from EvaluatorSecondary
void
AdditiveEvaluator::Evaluate_(const State& S,
                             CompositeVector& result)
{
  result.PutScalar(0.);
  
  for (std::map<Key, double>::const_iterator it=coefs_.begin();
       it!=coefs_.end(); ++it) {
    const CompositeVector& dep = S.Get<CompositeVector>(it->first, my_tag_);
    result.Update(it->second, dep, 1.0);
  }
}

void
AdditiveEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key,  const Key& wrt_tag,  CompositeVector& result)
{
  ASSERT(my_tag_==wrt_tag);
  result.PutScalar(coefs_[wrt_key]);
}


} // namespace
} // namespace

