/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The thaw depth evaluator gets the subsurface temperature.
  This computes the thaw depth over time.
  This is EvaluatorSecondaryMonotypeCV and depends on the subsurface temperature, 

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_MOISTURECONTENT_EVALUATOR_
#define AMANZI_FLOWRELATIONS_MOISTURECONTENT_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class MoistureContentEvaluator : public EvaluatorSecondaryMonotypeCV {

public:
  explicit
  MoistureContentEvaluator(Teuchos::ParameterList& plist);
  MoistureContentEvaluator(const MoistureContentEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;
  
protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);
  
    
  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);
  
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);


  bool updated_once_;
  bool volumetric_wc_,average_sat_;
  Key temp_key_, cv_key_, sat_key_, por_key_;
  Key domain_;
  double trans_width_;
private:
  static Utils::RegisteredFactory<Evaluator,MoistureContentEvaluator> reg_;

};
  
} //namespace
} //namespace 

#endif
