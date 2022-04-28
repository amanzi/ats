/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Author: Daniil Svyatsky (dasvyat@lanl.gov)
*/

//!
//! 
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

#ifndef AMANZI_FLOW_RELATIONS_DISTTILES_EVALUATOR_HH_
#define AMANZI_FLOW_RELATIONS_DISTTILES_EVALUATOR_HH_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class DistributedTilesRateEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  DistributedTilesRateEvaluator(Teuchos::ParameterList& plist);
  DistributedTilesRateEvaluator(const DistributedTilesRateEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override {
    return Teuchos::rcp(new DistributedTilesRateEvaluator(*this));
  }

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;
  

 protected:

  void InitializeFromPlist_();
  
  std::vector<double> times_;
  
  Key subsurface_marks_key_, sources_key_, pres_key_;
  Key domain_, domain_surf_;
  bool compatibility_checked_, implicit_;
  double p_enter_, k_;
  int num_ditches_;
  
 private:
  static Utils::RegisteredFactory<FieldEvaluator,DistributedTilesRateEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace
#endif
