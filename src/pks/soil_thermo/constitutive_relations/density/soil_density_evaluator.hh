/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

Evaluator for soil density.
----------------------------------------------------------------------------- */


#ifndef AMANZI_SOIL_DENSITY_EVALUATOR_HH_
#define AMANZI_SOIL_DENSITY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SoilThermo {

class SoilDensityEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  SoilDensityEvaluator(Teuchos::ParameterList& plist);
  SoilDensityEvaluator(const SoilDensityEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from SecondaryVariableEvaluator
  virtual void Evaluate_(const State& S,
      const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
      const Key& wrt_key, const Tag& wrt_tag,
      const std::vector<CompositeVector*>& result) override;

 protected:

  Key temperature_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,SoilDensityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
