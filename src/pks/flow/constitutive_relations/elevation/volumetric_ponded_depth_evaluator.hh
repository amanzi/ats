/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

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

#pragma once

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class VolumetricPondedDepthEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit VolumetricPondedDepthEvaluator(Teuchos::ParameterList& plist);
  VolumetricPondedDepthEvaluator(const VolumetricPondedDepthEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new VolumetricPondedDepthEvaluator(*this));
  }

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  Key pd_key_;
  Key delta_max_key_;
  Key delta_ex_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, VolumetricPondedDepthEvaluator> reg_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
