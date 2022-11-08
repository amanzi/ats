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

// TPLs
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"


#include "Factory.hh"
#include "Epetra_Vector_Factory.hh"
#include "EvaluatorSecondaryMonotype.hh" 

namespace Amanzi {
namespace Flow {
namespace Relations {

class DistributedTilesRateEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  DistributedTilesRateEvaluator(Teuchos::ParameterList& plist);
  DistributedTilesRateEvaluator(const DistributedTilesRateEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new DistributedTilesRateEvaluator(*this));
  }

  // Required methods from EvaluatorSecondaryMonotype
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>&  result) override;
  virtual void EnsureCompatibility_Structure_(State& S) override;
  

 protected:

  void InitializeFromPlist_();
  
  std::vector<double> times_;
  
  Key subsurface_marks_key_, dist_sources_key_, pres_key_, mol_dens_key_, factor_key_;
  Key domain_, domain_surf_;
  bool compatibility_checked_, implicit_;
  double p_enter_, k_;
  int num_ditches_, num_component_;
  Teuchos::RCP<Epetra_Vector> dist_src_vec_;
  
  
 private:
  static Utils::RegisteredFactory<Evaluator,DistributedTilesRateEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace
#endif
