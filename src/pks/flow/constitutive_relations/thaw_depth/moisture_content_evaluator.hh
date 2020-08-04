/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The thaw depth evaluator gets the subsurface temperature.
  This computes the thaw depth over time.
  This is SecondaryVariablesFieldEvaluator and depends on the subsurface temperature, 

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_MOISTURECONTENT_EVALUATOR_
#define AMANZI_FLOWRELATIONS_MOISTURECONTENT_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class MoistureContentEvaluator : public SecondaryVariableFieldEvaluator {

public:
  explicit
  MoistureContentEvaluator(Teuchos::ParameterList& plist);
  MoistureContentEvaluator(const MoistureContentEvaluator& other);
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
  bool volumetric_wc_;
  Key temp_key_, cv_key_, sat_key_, por_key_;
  Key domain_;
  double trans_width_;
private:
  static Utils::RegisteredFactory<FieldEvaluator,MoistureContentEvaluator> reg_;

};
  
} //namespace
} //namespace 

#endif
