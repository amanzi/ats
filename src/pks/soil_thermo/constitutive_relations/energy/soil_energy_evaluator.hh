/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

Evaluator for energy, e = cp*rho*T.
----------------------------------------------------------------------------- */


#ifndef AMANZI_SOIL_ENERGY_EVALUATOR_HH_
#define AMANZI_SOIL_ENERGY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SoilThermo {

class SoilEnergyEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  SoilEnergyEvaluator(Teuchos::ParameterList& plist);
  SoilEnergyEvaluator(const SoilEnergyEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from SecondaryVariableEvaluator
  virtual void Evaluate_(const State& S,
      const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
      const Key& wrt_key, const Tag& wrt_tag,
      const std::vector<CompositeVector*>& result) override;

 protected:

  Key temperature_key_;
  Key pres_key_;
  Key density_key_;
  Key heat_capacity_key_;
  Key pressure_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,SoilEnergyEvaluator> factory_;

};

} // namespace
} // namespace

#endif
