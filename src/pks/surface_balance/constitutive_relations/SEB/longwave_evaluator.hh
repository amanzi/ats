/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Author: Ethan Coon (ecoon@ornl.gov)
*/

//! Evaluates incoming longwave radiation from rel humidity and air temperature.

/*!

Requires the following dependencies:

* `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
* `"relative humicity key`" ``[string]`` **DOMAIN-relative_humidity**
         
*/

#ifndef AMANZI_SURFACE_BALANCE_LONGWAVE_EVALUATOR_HH_
#define AMANZI_SURFACE_BALANCE_LONGWAVE_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {

class LongwaveEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  LongwaveEvaluator(Teuchos::ParameterList& plist);
  LongwaveEvaluator(const LongwaveEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new LongwaveEvaluator(*this));
  }

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Exceptions::amanzi_throw("NotImplemented: LongwaveEvaluator currently does not provide derivatives.");
  }

 protected:

  Key air_temp_key_, rel_hum_key_;
  double min_rel_hum_;
  double stephB_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LongwaveEvaluator> reg_;

};

} //namespace
} //namespace

#endif
