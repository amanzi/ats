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
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class SurfDistributedTilesRateEvaluator : public EvaluatorSecondary {

 public:
  explicit
  SurfDistributedTilesRateEvaluator(Teuchos::ParameterList& plist);
  SurfDistributedTilesRateEvaluator(const SurfDistributedTilesRateEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new SurfDistributedTilesRateEvaluator(*this));
  }

  virtual void EnsureCompatibility(State& S) override;

 protected:

  virtual void Update_(State& S) override ;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag)  override {};
  
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
