/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The thaw depth evaluator gets the subsurface temperature.
  This computes the thaw depth over time.
  This is SecondaryVariablesFieldEvaluator and depends on the subsurface temperature, 

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_THAWDEPTH_EVALUATOR_
#define AMANZI_FLOWRELATIONS_THAWDEPTH_EVALUATOR_

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Flow {

class ThawDepthEvaluator : public SecondaryVariableFieldEvaluator {

public:
  explicit
  ThawDepthEvaluator(Teuchos::ParameterList& plist);
  ThawDepthEvaluator(const ThawDepthEvaluator& other);
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
  double trans_width_;
  Key domain_;
  Key temp_key_;
private:
  static Utils::RegisteredFactory<FieldEvaluator,ThawDepthEvaluator> reg_;

};
  
} //namespace
} //namespace 

#endif
