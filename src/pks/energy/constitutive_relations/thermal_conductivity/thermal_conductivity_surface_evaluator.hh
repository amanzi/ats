/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_RELATIONS_TC_SURFACE_EVALUATOR_HH_
#define AMANZI_ENERGY_RELATIONS_TC_SURFACE_EVALUATOR_HH_

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Energy {

class ThermalConductivitySurfaceEvaluator :
    public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  ThermalConductivitySurfaceEvaluator(Teuchos::ParameterList& plist);
  ThermalConductivitySurfaceEvaluator(const ThermalConductivitySurfaceEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from SecondaryVariableFieldModel
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

 protected:
  // dependencies
  Key uf_key_;
  Key height_key_;

  double K_liq_;
  double K_ice_;
  double min_K_;
};

} // namespace
} // namespace

#endif
