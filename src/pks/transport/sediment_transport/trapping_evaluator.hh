/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The erosion evaluator gets the erosion rates.


  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_TRAPPINGRATE_EVALUATOR_
#define AMANZI_TRAPPINGRATE_EVALUATOR_

// TPLs
#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {

class TrappingRateEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  TrappingRateEvaluator(Teuchos::ParameterList& plist);

  TrappingRateEvaluator(const TrappingRateEvaluator& other);
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

  double visc_, d_p_, alpha_, beta_, gamma_;

  Key velocity_key_;
  Key sediment_key_;
  Key ponded_depth_key_;
  Key biomass_key_;  
  double sediment_density_;
  static Utils::RegisteredFactory<FieldEvaluator,TrappingRateEvaluator> factory_;

};

} //namespace

#endif
