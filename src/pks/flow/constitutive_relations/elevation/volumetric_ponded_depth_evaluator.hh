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

#pragma once

#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class VolumetricPondedDepthEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  VolumetricPondedDepthEvaluator(Teuchos::ParameterList& plist);
  VolumetricPondedDepthEvaluator(const VolumetricPondedDepthEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new VolumetricPondedDepthEvaluator(*this));
  }

 protected:
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key pd_key_;
  Key delta_max_key_;
  Key delta_ex_key_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,VolumetricPondedDepthEvaluator> reg_;

};

} //namespace
} //namespace


