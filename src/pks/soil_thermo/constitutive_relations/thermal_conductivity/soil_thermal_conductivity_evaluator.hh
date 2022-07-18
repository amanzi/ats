/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SOIL_TC_EVALUATOR_HH_
#define AMANZI_SOIL_TC_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SoilThermo {

class SoilThermalConductivityEvaluator :
    public EvaluatorSecondaryMonotypeCV {

 public:
  // constructor format for all derived classes
  SoilThermalConductivityEvaluator(Teuchos::ParameterList& plist);
  SoilThermalConductivityEvaluator(const SoilThermalConductivityEvaluator& other);

  Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from SecondaryVariableFieldModel
  virtual void Evaluate_(const State& S,
      const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
      const Key& wrt_key, const Tag& wrt_tag,
      const std::vector<CompositeVector*>& result) override;

 protected:
  // dependencies

  Key temperature_key_;
  Key water_content_key_;
  Key ice_content_key_;
  Key cell_is_ice_key_;
  Key cell_vol_key_;
  Key pressure_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,SoilThermalConductivityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
