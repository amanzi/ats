/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The gas content evaluator gets the subsurface temperature.
  This computes the gas content in the underneath column and puts on the surface cell.
  This is SecondaryVariablesFieldEvaluator and depends on the subsurface temperature, 

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_GASCONTENT_EVALUATOR_
#define AMANZI_FLOWRELATIONS_GASCONTENT_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class GasContentEvaluator : public SecondaryVariableFieldEvaluator {

public:
  explicit
  GasContentEvaluator(Teuchos::ParameterList& plist);
  GasContentEvaluator(const GasContentEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;
  
protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);
  
    
  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);
  
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);


  bool updated_once_;
  Key temp_key_, cv_key_, sat_key_;
  Key domain_;
  double trans_width_;
private:
  static Utils::RegisteredFactory<FieldEvaluator,GasContentEvaluator> reg_;

};
  
} //namespace
} //namespace 

#endif
