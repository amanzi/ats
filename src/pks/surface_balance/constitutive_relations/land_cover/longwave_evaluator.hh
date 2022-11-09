/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Evaluates incoming longwave radiation from vapor pressure and air temperature.

/*!

.. _longwave_evaluator-spec:
.. admonition:: longwave_evaluator-spec


    DEPENDENCIES:

    * `"air temperature key`" ``[string]`` **DOMAIN-air_temperature**
    * `"vapor pressure air key`" ``[string]`` **DOMAIN-vapor_pressure_air**

*/

#ifndef AMANZI_SURFACE_BALANCE_LONGWAVE_EVALUATOR_HH_
#define AMANZI_SURFACE_BALANCE_LONGWAVE_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class LongwaveEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit LongwaveEvaluator(Teuchos::ParameterList& plist);
  LongwaveEvaluator(const LongwaveEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new LongwaveEvaluator(*this));
  }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override {
    Exceptions::amanzi_throw("NotImplemented: LongwaveEvaluator currently does not provide derivatives.");
  }

 protected:
  Key air_temp_key_, vp_air_key_;
  double scale_;

 private:
  static Utils::RegisteredFactory<Evaluator,LongwaveEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
