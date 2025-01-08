/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  This WRM model evaluates the saturation of ice, water, and gas.

*/

#ifndef AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_WRM_PERMAFROST_EVALUATOR_

#include "wrm.hh"
#include "wrm_partition.hh"
#include "wrm_permafrost_model.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMPermafrostEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit WRMPermafrostEvaluator(Teuchos::ParameterList& plist);
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrms);
  WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<WRMPermafrostModelPartition>& models);
  WRMPermafrostEvaluator(const WRMPermafrostEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }
  Teuchos::RCP<WRMPermafrostModelPartition> get_WRMPermafrostModels() { return permafrost_models_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void EnsureCompatibility_Structure_(State& S) override
  {
    EnsureCompatibility_StructureSame_(S);
  }

  void InitializeFromPlist_();

 protected:
  Key pc_liq_key_;
  Key pc_ice_key_;

  Teuchos::RCP<WRMPermafrostModelPartition> permafrost_models_;
  Teuchos::RCP<WRMPartition> wrms_;

 private:
  static Utils::RegisteredFactory<Evaluator, WRMPermafrostEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
