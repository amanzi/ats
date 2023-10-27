/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Capillary pressure for gas on a liquid.
/*!

.. _pc-liquid-evaluator-spec:
.. admonition:: pc-liquid-evaluator-spec

   KEYS:

   - `"pressure`" **DOMAIN-pressure**

*/


#ifndef AMANZI_RELATIONS_PC_LIQUID_EVALUATOR_HH_
#define AMANZI_RELATIONS_PC_LIQUID_EVALUATOR_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class PCLiqAtm;

class PCLiquidEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit PCLiquidEvaluator(Teuchos::ParameterList& plist);
  PCLiquidEvaluator(const PCLiquidEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<PCLiqAtm> get_PCLiqAtm() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  // the actual model
  Teuchos::RCP<PCLiqAtm> model_;

  // Keys for fields
  // dependencies
  Key pres_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, PCLiquidEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
