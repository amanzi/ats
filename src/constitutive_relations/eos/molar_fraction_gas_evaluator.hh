/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Determining the molar fraction of a gas component within a gas mixture.

*/

#ifndef AMANZI_RELATIONSRELATIONS_MOLAR_FRACTION_GAS_
#define AMANZI_RELATIONSRELATIONS_MOLAR_FRACTION_GAS_

#include "vapor_pressure_relation.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class MolarFractionGasEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit MolarFractionGasEvaluator(Teuchos::ParameterList& plist);
  MolarFractionGasEvaluator(const MolarFractionGasEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  Teuchos::RCP<VaporPressureRelation> get_VaporPressureRelation() { return sat_vapor_model_; }

 protected:
  Key temp_key_;

  Teuchos::RCP<VaporPressureRelation> sat_vapor_model_;

 private:
  static Utils::RegisteredFactory<Evaluator, MolarFractionGasEvaluator> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
