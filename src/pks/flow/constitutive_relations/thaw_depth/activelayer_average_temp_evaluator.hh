/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The Active layer average temperature evaluator gets the subsurface temperature.
  This computes the average active layer temperature.
  This is EvaluatorSecondaryMonotypeCV and depends on the subsurface temperature,

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_ALTTEMP_EVALUATOR_
#define AMANZI_FLOWRELATIONS_ALTTEMP_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ActiveLayerAverageTempEvaluator : public EvaluatorSecondaryMonotypeCV {

public:
  explicit
  ActiveLayerAverageTempEvaluator(Teuchos::ParameterList& plist);
  ActiveLayerAverageTempEvaluator(const ActiveLayerAverageTempEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  virtual bool Update(State& S, const Key& request) override;
  virtual void EnsureCompatibility(State& S) override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
                         const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
               const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override;

 protected:
  bool updated_once_;
  Key temp_key_;
  Key domain_;
  double trans_width_;

private:
  static Utils::RegisteredFactory<Evaluator,ActiveLayerAverageTempEvaluator> reg_;

};

} //namespace
} //namespace

#endif
