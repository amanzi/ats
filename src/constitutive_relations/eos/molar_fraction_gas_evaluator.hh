/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONSRELATIONS_MOLAR_FRACTION_GAS_
#define AMANZI_RELATIONSRELATIONS_MOLAR_FRACTION_GAS_

#include "vapor_pressure_relation.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class MolarFractionGasEvaluator : public EvaluatorSecondary<CompositeVector, CompositeVectorSpace> {

 public:
  explicit
  MolarFractionGasEvaluator(Teuchos::ParameterList& plist);

  MolarFractionGasEvaluator(const MolarFractionGasEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondary
  virtual void Evaluate_(const State& S,
                         CompositeVector& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key, const Key& wrt_tag,
                                          CompositeVector& result) override;

  Teuchos::RCP<VaporPressureRelation> get_VaporPressureRelation() {
    return sat_vapor_model_; }

 protected:
  Key temp_key_;

  Teuchos::RCP<VaporPressureRelation> sat_vapor_model_;

 private:
  static Utils::RegisteredFactory<Evaluator,MolarFractionGasEvaluator> factory_;

};

} //namespace
} //namespace

#endif
