/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Determine the volumetric ponded depth from ponded depth and subgrid parameters.
/*!

* `"maximum relief key`" ``[string]`` **DOMAIN-maximum_relief**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"ponded depth key`" ``[string]`` **DOMAIN-ponded_depth**
         The true height of the water surface.

*/

#include "volumetric_ponded_depth_evaluator.hh"
#include "subgrid_microtopography.hh"

namespace Amanzi {
namespace Flow {

VolumetricPondedDepthEvaluator::VolumetricPondedDepthEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // dependencies
  pd_key_ = Keys::readKey(plist_, domain, "ponded depth key", "ponded_depth");
  dependencies_.insert(KeyTag{pd_key_, tag});

  delta_max_key_ = Keys::readKey(plist_, domain, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(KeyTag{delta_max_key_, tag});

  delta_ex_key_ = Keys::readKey(plist_, domain, "excluded volume", "excluded_volume");
  dependencies_.insert(KeyTag{delta_ex_key_, tag});
}


void
VolumetricPondedDepthEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  for (const auto& comp : *result[0]) {
    auto& res = *result[0]->ViewComponent(comp,false);
    const auto& pd = *S.Get<CompositeVector>(pd_key_, tag).ViewComponent(comp,false);
    const auto& del_max = *S.Get<CompositeVector>(delta_max_key_, tag).ViewComponent(comp,false);
    const auto& del_ex = *S.Get<CompositeVector>(delta_ex_key_, tag).ViewComponent(comp,false);

    for (int c=0; c!=res.MyLength(); ++c){
      AMANZI_ASSERT(Microtopography::validParameters(del_max[0][c], del_ex[0][c]));
      res[0][c] = Microtopography::volumetricDepth(pd[0][c], del_max[0][c], del_ex[0][c]);
    }
  }
}


void
VolumetricPondedDepthEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  for (const auto& comp : *result[0]) {
    auto& res = *result[0]->ViewComponent(comp,false);
    const auto& pd = *S.Get<CompositeVector>(pd_key_, tag).ViewComponent(comp,false);
    const auto& del_max = *S.Get<CompositeVector>(delta_max_key_, tag).ViewComponent(comp,false);
    const auto& del_ex = *S.Get<CompositeVector>(delta_ex_key_, tag).ViewComponent(comp,false);

    if (wrt_key == pd_key_) {
      for (int c=0; c!=res.MyLength(); ++c){
        res[0][c] = Microtopography::dVolumetricDepth_dDepth(pd[0][c], del_max[0][c], del_ex[0][c]);
        res[0][c] = std::max(res[0][c],0.001);
      }
    } else {
      Errors::Message msg("VolumetricPondedDepthEvaluator: Not Implemented: no derivatives implemented other than ponded depth.");
      Exceptions::amanzi_throw(msg);
    }
  }
}


} //namespace
} //namespace
