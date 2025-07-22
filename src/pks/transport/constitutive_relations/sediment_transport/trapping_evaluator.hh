/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  The erosion evaluator gets the erosion rates.


*/
#pragma once

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class TrappingRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit TrappingRateEvaluator(Teuchos::ParameterList& plist);

  TrappingRateEvaluator(const TrappingRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  double visc_, d_p_, alpha_, beta_, gamma_;

  Key velocity_key_;
  Key sediment_key_;
  Key ponded_depth_key_;
  Key stem_diameter_key_, stem_height_key_, stem_density_key_;
  double sediment_density_;

  static Utils::RegisteredFactory<Evaluator, TrappingRateEvaluator> factory_;
};

} // namespace Amanzi
