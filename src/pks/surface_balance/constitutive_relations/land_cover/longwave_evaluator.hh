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
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

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

  Key air_temp_key_, vp_air_key_;
  double scale_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,LongwaveEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
