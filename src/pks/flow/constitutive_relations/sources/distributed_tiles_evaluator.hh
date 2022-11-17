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

#ifndef AMANZI_FLOW_RELATIONS_DISTTILES_EVALUATOR_HH_
#define AMANZI_FLOW_RELATIONS_DISTTILES_EVALUATOR_HH_

// TPLs
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"


#include "Factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class DistributedTilesRateEvaluator : public EvaluatorSecondary {

 public:
  explicit
  DistributedTilesRateEvaluator(Teuchos::ParameterList& plist);
  DistributedTilesRateEvaluator(const DistributedTilesRateEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new DistributedTilesRateEvaluator(*this));
  }

  virtual void EnsureCompatibility(State& S) override;
  
 protected:
  // Required methods from EvaluatorSecondary
  virtual void Update_(State& S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag)  override {};


  Key domain_, domain_surf_;

  Key subsurface_marks_key_;
  Key dist_sources_key_;
  Key pres_key_;
  Key mol_dens_key_;
  Key factor_key_;

  bool compatibility_checked_, implicit_;

  double p_enter_;
  double k_;
  int num_ditches_;
  int num_component_;

 private:
  static Utils::RegisteredFactory<Evaluator,DistributedTilesRateEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace
#endif
