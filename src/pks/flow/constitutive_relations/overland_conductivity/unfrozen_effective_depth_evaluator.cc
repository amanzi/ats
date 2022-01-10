/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen effective depth.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_effective_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

UnfrozenEffectiveDepthEvaluator::UnfrozenEffectiveDepthEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  Key domain = Keys::getDomain(my_key_);

  depth_key_ = plist_.get<std::string>("depth key", Keys::getKey(domain,"ponded_depth"));
  dependencies_.insert(depth_key_);

  uf_key_ = plist_.get<std::string>("unfrozen fraction key", Keys::getKey(domain,"unfrozen_fraction"));
  alpha_ = plist_.get<double>("ice retardation exponent [-]", 1.0);
  dependencies_.insert(uf_key_);

  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("unfrozen effective depth key",
            "unfrozen_effective_depth");
  }

}


Teuchos::RCP<FieldEvaluator>
UnfrozenEffectiveDepthEvaluator::Clone() const {
  return Teuchos::rcp(new UnfrozenEffectiveDepthEvaluator(*this));
}


// Required methods from SecondaryVariableFieldEvaluator
void UnfrozenEffectiveDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> depth = S->GetPtr<CompositeVector>(depth_key_);
  Teuchos::RCP<const CompositeVector> uf = S->GetPtr<CompositeVector>(uf_key_);
  
  for (auto compname : *result) {
    auto& result_c = *result->ViewComponent(compname, false);
    const auto& depth_c = *depth->ViewComponent(compname, false);
    const auto& uf_c = *uf->ViewComponent(compname, false);

    for (int c=0; c!=result_c.MyLength(); ++c) {
      result_c[0][c] = depth_c[0][c] * std::pow(uf_c[0][c], alpha_);
    }
  }
}


void UnfrozenEffectiveDepthEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S,
    Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> depth = S->GetPtr<CompositeVector>(depth_key_);
  Teuchos::RCP<const CompositeVector> uf = S->GetPtr<CompositeVector>(uf_key_);

  if (wrt_key == depth_key_) {
    for (auto compname : *result) {
      auto& result_c = *result->ViewComponent(compname, false);
      const auto& depth_c = *depth->ViewComponent(compname, false);
      const auto& uf_c = *uf->ViewComponent(compname, false);

      for (int c=0; c!=result_c.MyLength(); ++c) {
        result_c[0][c] = std::pow(uf_c[0][c], alpha_);
      }
    }
  } else if (wrt_key == uf_key_) {
    for (auto compname : *result) {
      auto& result_c = *result->ViewComponent(compname, false);
      const auto& depth_c = *depth->ViewComponent(compname, false);
      const auto& uf_c = *uf->ViewComponent(compname, false);

      for (int c=0; c!=result_c.MyLength(); ++c) {
        result_c[0][c] = alpha_ * depth_c[0][c] * std::pow(uf_c[0][c], alpha_ - 1);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }

  
}


} //namespace
} //namespace

