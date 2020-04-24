/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The erosion evaluator gets the erosion rates.


  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_SETTLEMENTRATE_EVALUATOR_
#define AMANZI_SETTLEMENTRATE_EVALUATOR_

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {

class SettlementRateEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  SettlementRateEvaluator(Teuchos::ParameterList& plist);

  SettlementRateEvaluator(const SettlementRateEvaluator& other);
  virtual Teuchos::RCP<FieldEvaluator> Clone() const;
  
  // virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
  //         const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;

  // virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);

  //virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S){};

protected:

    // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                               Key wrt_key,
                                               const Teuchos::Ptr<CompositeVector>& result);

  double tau_d_;
  double ws_;
  double gamma_;
  double lambda_, umax_, xi_;
  double sediment_density_;
  Key velocity_key_, sediment_key_;

  static Utils::RegisteredFactory<FieldEvaluator,SettlementRateEvaluator> factory_;

};

} //namespace

#endif
