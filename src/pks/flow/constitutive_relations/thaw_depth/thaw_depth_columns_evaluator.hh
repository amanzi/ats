/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The thaw depth evaluator gets the subsurface temperature.
  This computes the thaw depth over time.
  This is EvaluatorSecondaryMonotypeCV and depends on the subsurface temperature, 

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_THAWDEPTH_EVALUATOR_
#define AMANZI_FLOWRELATIONS_THAWDEPTH_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ThawDepthColumnsEvaluator : public EvaluatorSecondaryMonotypeCV {

public:
  explicit
  ThawDepthColumnsEvaluator(Teuchos::ParameterList& plist);
  ThawDepthColumnsEvaluator(const ThawDepthColumnsEvaluator& other);
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
  double trans_width_;
  Key domain_;
  Key temp_key_;
private:
  static Utils::RegisteredFactory<Evaluator,ThawDepthColumnsEvaluator> reg_;

};
  
} //namespace
} //namespace 

#endif
