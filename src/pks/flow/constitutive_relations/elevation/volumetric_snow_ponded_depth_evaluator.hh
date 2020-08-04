/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Determine the volumetric ponded and snow depths from ponded depth and snow depth.
/*!

* `"maximum relief key`" ``[string]`` **DOMAIN-maximum_relief**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"ponded depth key`" ``[string]`` **DOMAIN-ponded_depth**
         The true height of the water surface.
* `"snow depth key`" ``[string]`` **DOMAIN-snow_depth**
         The true height of the water surface.

*/

#pragma once

#include "secondary_variables_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {


class VolumetricSnowPondedDepthEvaluator : public SecondaryVariablesFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  VolumetricSnowPondedDepthEvaluator(Teuchos::ParameterList& plist);
  VolumetricSnowPondedDepthEvaluator(const VolumetricSnowPondedDepthEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override {
    return Teuchos::rcp(new VolumetricSnowPondedDepthEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 protected:
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;

  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results) override;

 protected:

  Key vol_pd_key_;
  Key vol_sd_key_;
  Key pd_key_;
  Key sd_key_;
  Key delta_max_key_;
  Key delta_ex_key_;

  Key domain_snow_, domain_surf_;

  bool compatibility_checked_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,VolumetricSnowPondedDepthEvaluator> reg_;

};

} //namespace
} //namespace


