/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  Suction head = \Psi( sat )

*/

#ifndef AMANZI_FLOWRELATIONS_SUCTION_HEAD_EVALUATOR_
#define AMANZI_FLOWRELATIONS_SUCTION_HEAD_EVALUATOR_

#include "wrm.hh"
#include "wrm_partition.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {            
namespace Flow {

class SuctionHeadEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit SuctionHeadEvaluator(Teuchos::ParameterList& plist);
  SuctionHeadEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrms);
  SuctionHeadEvaluator(const SuctionHeadEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  void InitializeFromPlist_();

  Teuchos::RCP<WRMPartition> wrms_;
  Key sat_key_;
  double min_val_;

 private:
  static Utils::RegisteredFactory<Evaluator, SuctionHeadEvaluator> factory_;
};

} // namespace Flow
}  
} // namespace Amanzi

#endif
