/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The erosion evaluator gets the erosion rates.


  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#pragma once

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class OrganicMatterRateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit OrganicMatterRateEvaluator(Teuchos::ParameterList& plist);
  OrganicMatterRateEvaluator(const OrganicMatterRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;


  double Bmax_;
  double Q_db0_;
  double Q_on_Bmax_;
  Key biomass_key_;

  static Utils::RegisteredFactory<Evaluator, OrganicMatterRateEvaluator> factory_;
};

} // namespace Amanzi
