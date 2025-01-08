/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Evaluates the unfrozen effective depth.

*/

#include "unfrozen_effective_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

UnfrozenEffectiveDepthEvaluator::UnfrozenEffectiveDepthEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(my_keys_.front().first);

  depth_key_ = Keys::readKey(plist_, domain, "depth", "ponded_depth");
  dependencies_.insert(KeyTag{ depth_key_, tag });

  uf_key_ = Keys::readKey(plist_, domain, "unfrozen fraction", "unfrozen_fraction");
  dependencies_.insert(KeyTag{ uf_key_, tag });

  alpha_ = plist_.get<double>("ice retardation exponent [-]", 1.0);
}


Teuchos::RCP<Evaluator>
UnfrozenEffectiveDepthEvaluator::Clone() const
{
  return Teuchos::rcp(new UnfrozenEffectiveDepthEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
UnfrozenEffectiveDepthEvaluator::Evaluate_(const State& S,
                                           const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> depth = S.GetPtr<CompositeVector>(depth_key_, tag);
  Teuchos::RCP<const CompositeVector> uf = S.GetPtr<CompositeVector>(uf_key_, tag);

  for (auto compname : *result[0]) {
    auto& result_c = *result[0]->ViewComponent(compname, false);
    const auto& depth_c = *depth->ViewComponent(compname, false);
    const auto& uf_c = *uf->ViewComponent(compname, false);

    for (int c = 0; c != result_c.MyLength(); ++c) {
      result_c[0][c] = depth_c[0][c] * std::pow(uf_c[0][c], alpha_);
    }
  }
}


void
UnfrozenEffectiveDepthEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> depth = S.GetPtr<CompositeVector>(depth_key_, tag);
  Teuchos::RCP<const CompositeVector> uf = S.GetPtr<CompositeVector>(uf_key_, tag);

  if (wrt_key == depth_key_) {
    for (auto compname : *result[0]) {
      auto& result_c = *result[0]->ViewComponent(compname, false);
      const auto& depth_c = *depth->ViewComponent(compname, false);
      const auto& uf_c = *uf->ViewComponent(compname, false);

      for (int c = 0; c != result_c.MyLength(); ++c) {
        result_c[0][c] = std::pow(uf_c[0][c], alpha_);
      }
    }
  } else if (wrt_key == uf_key_) {
    for (auto compname : *result[0]) {
      auto& result_c = *result[0]->ViewComponent(compname, false);
      const auto& depth_c = *depth->ViewComponent(compname, false);
      const auto& uf_c = *uf->ViewComponent(compname, false);

      for (int c = 0; c != result_c.MyLength(); ++c) {
        result_c[0][c] = alpha_ * depth_c[0][c] * std::pow(uf_c[0][c], alpha_ - 1);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}


} // namespace Flow
} // namespace Amanzi
