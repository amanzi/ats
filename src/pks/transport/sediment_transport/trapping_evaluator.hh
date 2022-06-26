/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The erosion evaluator gets the erosion rates.


  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_TRAPPINGRATE_EVALUATOR_
#define AMANZI_TRAPPINGRATE_EVALUATOR_

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {

class TrappingRateEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  TrappingRateEvaluator(Teuchos::ParameterList& plist);

  TrappingRateEvaluator(const TrappingRateEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;
  
  // virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
  //         const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;

  // virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);

  //virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S){};

  protected:


    // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override;  

  double visc_, d_p_, alpha_, beta_, gamma_;

  Key velocity_key_;
  Key sediment_key_;
  Key ponded_depth_key_;
  Key biomass_key_;  
  double sediment_density_;
  static Utils::RegisteredFactory<Evaluator,TrappingRateEvaluator> factory_;

};

} //namespace

#endif
