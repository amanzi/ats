/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The column average temperature evaluator gets the subsurface temperature.
  This computes the average column temperature to a specified depth.
  This is EvaluatorSecondaryMonotypeCV and depends on the subsurface temperature,

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_COLUMNTEMP_EVALUATOR_
#define AMANZI_FLOWRELATIONS_COLUMNTEMP_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ColumnAverageTempEvaluator : public EvaluatorSecondaryMonotypeCV {

public:
  explicit
  ColumnAverageTempEvaluator(Teuchos::ParameterList& plist);
  ColumnAverageTempEvaluator(const ColumnAverageTempEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;


  virtual bool Update(State& S, const Key& request) override;
  virtual void EnsureCompatibility(State& S) override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
                         const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& result) override;

 protected:
  bool updated_once_;
  Key temp_key_;
  Key domain_;
  int ncells_depth_;
  double depth_;
private:
  static Utils::RegisteredFactory<Evaluator,ColumnAverageTempEvaluator> reg_;

};

} //namespace
} //namespace

#endif
