/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

Thermal conductivity based on two-phases (air,liquid) composition of the
porous media.

`"evaluator type`" = `"two-phase thermal conductivity`"

.. _evaluator-two-phase-thermal-conductivity-spec:
.. admonition:: evaluator-two-phase-thermal-conductivity-spec

   * `"thermal conductivity parameters`" ``[thermal-conductivity-twophase-typed-spec-list]``

   KEYS:

   - `"porosity`"
   - `"saturation liquid`"

*/

#ifndef AMANZI_ENERGY_RELATIONS_TC_TWOPHASE_EVALUATOR_HH_
#define AMANZI_ENERGY_RELATIONS_TC_TWOPHASE_EVALUATOR_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "thermal_conductivity_twophase.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

// Equation of State model
class ThermalConductivityTwoPhaseEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  using RegionModelPair = std::pair<std::string, Teuchos::RCP<ThermalConductivityTwoPhase>>;

  // constructor format for all derived classes
  ThermalConductivityTwoPhaseEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivityTwoPhaseEvaluator(const ThermalConductivityTwoPhaseEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  std::vector<RegionModelPair> tcs_;

  // Keys for fields
  // dependencies
  Key poro_key_;
  Key sat_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, ThermalConductivityTwoPhaseEvaluator> factory_;
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
