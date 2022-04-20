/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_RELATIONS_TC_TWOPHASE_EVALUATOR_HH_
#define AMANZI_ENERGY_RELATIONS_TC_TWOPHASE_EVALUATOR_HH_

#include "secondary_variable_field_evaluator.hh"
#include "thermal_conductivity_twophase.hh"

namespace Amanzi {
namespace Energy {

// Equation of State model
class ThermalConductivityTwoPhaseEvaluator :
    public SecondaryVariableFieldEvaluator {

 public:

  typedef std::pair<std::string,Teuchos::RCP<ThermalConductivityTwoPhase> > RegionModelPair;

  // constructor format for all derived classes
  ThermalConductivityTwoPhaseEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivityTwoPhaseEvaluator(const ThermalConductivityTwoPhaseEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  // Required methods from SecondaryVariableFieldModel
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  std::vector<RegionModelPair> tcs_;

  // Keys for fields
  // dependencies
  Key poro_key_;
  Key sat_key_;
};

} // namespace
} // namespace

#endif
