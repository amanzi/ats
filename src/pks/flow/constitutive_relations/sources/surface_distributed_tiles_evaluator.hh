/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Author: Daniil Svyatsky (dasvyat@lanl.gov)
*/

//!
/*!

Requires the following dependencies:

*
*

Allows the following parameters:

*


*


*

.. note:

*/

#ifndef AMANZI_FLOW_RELATIONS_SURFDISTTILES_EVALUATOR_HH_
#define AMANZI_FLOW_RELATIONS_SURFDISTTILES_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class SurfDistributedTilesRateEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  SurfDistributedTilesRateEvaluator(Teuchos::ParameterList& plist);
  SurfDistributedTilesRateEvaluator(const SurfDistributedTilesRateEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new SurfDistributedTilesRateEvaluator(*this));
  }

 protected:
  // Required methods from EvaluatorSecondaryMonotype
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>&  result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  std::vector<double> times_;
  Key surface_marks_key_, surf_len_key_, dist_sources_key_;
  Key domain_, domain_surf_;
  bool compatibility_checked_, implicit_;
  int num_ditches_;

 private:
  static Utils::RegisteredFactory<Evaluator,SurfDistributedTilesRateEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace
#endif
