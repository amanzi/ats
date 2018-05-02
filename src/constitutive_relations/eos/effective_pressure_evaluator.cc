/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EffectivePressureEvaluator evaluates p_eff = max(p_atm, p_liquid), which is used for EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "factory.hh"
#include "effective_pressure_evaluator.hh"

namespace Amanzi {
namespace Relations {

EffectivePressureEvaluator::EffectivePressureEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondary(plist) {
  if (my_key_ == std::string("")) {
    my_key_ = ep_plist_.get<std::string>("effective pressure key", "effective_pressure");
  }

  Key domain_name = Keys::getDomain(my_key_);

  // -- pressure
  pres_key_ = plist_.get<std::string>("pressure key",
          Keys::getKey(domain_name, "pressure"));
  dependencies_.emplace_back(std::make_pair(pres_key_, my_tag_));

  // -- logging
  if (vo_.os_OK(Teuchos::VERB_EXTREME)) {
    Teuchos::OSTab tab = vo_.getOSTab();
    for (auto& dep : dependencies_) {
      *vo_.os() << " dep: " << dep.first <<" "<<dep.second << std::endl;
    }
  }

}


EffectivePressureEvaluator::EffectivePressureEvaluator(
        const EffectivePressureEvaluator& other) :
    EvaluatorSecondary(other),
    pres_key_(other.pres_key_) {}


Teuchos::RCP<Evaluator> EffectivePressureEvaluator::Clone() const {
  return Teuchos::rcp(new EffectivePressureEvaluator(*this));
}

void EffectivePressureEvaluator::Evaluate_(const State& S,
                                           CompositeVector& result) {
  // Pull dependencies out of state.
  const CompositeVector& pres = S.Get<CompositeVector>(pres_key_, my_tag_);
  const double& p_atm = S.Get<double>("atmospheric_pressure");

  // evaluate effective pressure as max(pres, p_atm)
  for (CompositeVector::name_iterator comp=result.begin();
       comp!=result.end(); ++comp) {
    const Epetra_MultiVector& pres_v = *(pres.ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result.ViewComponent(*comp,false));

    int count = result.size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = std::max<double>(pres_v[0][id], p_atm);
    }
  }
}

void EffectivePressureEvaluator::EvaluatePartialDerivative_(const State& S,
                                                            const Key& wrt_key, const Key& wrt_tag, CompositeVector& result) {

   ASSERT(wrt_tag == my_tag_);
  
  // Pull dependencies out of state.
  const CompositeVector& pres = S.Get<CompositeVector>(pres_key_);
  const double& p_atm = S.Get<double>("atmospheric_pressure");

  ASSERT(wrt_key == pres_key_);
  // pressure is max(pres, p_atm), so derivative is 1 or 0
  for (CompositeVector::name_iterator comp=result.begin();
       comp!=result.end(); ++comp) {
    const Epetra_MultiVector& pres_v = *(pres.ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result.ViewComponent(*comp,false));

    int count = result.size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = pres_v[0][id] > p_atm ? 1.0 : 0.0;
    }
  }
};

} // namespace
} // namespace


