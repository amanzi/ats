/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The thaw depth evaluator gets the subsurface temperature.
  This computes the thaw depth over time.
  This is EvaluatorSecondaryMonotypeCV and depends on the subsurface temperature, 

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_WATERTABLE_EVALUATOR_
#define AMANZI_FLOWRELATIONS_WATERTABLE_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class WaterTableColumnsEvaluator : public EvaluatorSecondaryMonotypeCV {

public:
  explicit
  WaterTableColumnsEvaluator(Teuchos::ParameterList& plist);
  WaterTableColumnsEvaluator(const WaterTableColumnsEvaluator& other);
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
  Key temp_key_, sat_key_, pd_key_;
  Key domain_;
  double trans_width_;
private:
  static Utils::RegisteredFactory<Evaluator,WaterTableColumnsEvaluator> reg_;

};
  
} //namespace
} //namespace 

#endif
