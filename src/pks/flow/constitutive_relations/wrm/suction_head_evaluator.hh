/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Suction head = \Psi( sat )

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_SUCTION_HEAD_EVALUATOR_
#define AMANZI_FLOWRELATIONS_SUCTION_HEAD_EVALUATOR_

#include "wrm.hh"
#include "wrm_partition.hh"
#include "secondary_variable_field_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class SuctionHeadEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  SuctionHeadEvaluator(Teuchos::ParameterList& plist);

  SuctionHeadEvaluator(Teuchos::ParameterList& plist,
                   const Teuchos::RCP<WRMPartition>& wrms);

  SuctionHeadEvaluator(const SuctionHeadEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  
  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }

 protected:
  
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  void InitializeFromPlist_();

  Teuchos::RCP<WRMPartition> wrms_;
  Key sat_key_;  
  double min_val_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,SuctionHeadEvaluator> factory_;
};

} //namespace
} //namespace

#endif
