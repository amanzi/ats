/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The column average temperature evaluator gets the subsurface temperature.
  This computes the average column temperature to a specified depth.
  This is EvaluatorSecondaryMonotypeCV and depends on the subsurface temperature, 

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_COLUMNTEMP_EVALUATOR_
#define AMANZI_FLOWRELATIONS_COLUMNTEMP_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ColumnAverageTempEvaluator : public EvaluatorSecondaryMonotypeCV {

public:
  explicit
  ColumnAverageTempEvaluator(Teuchos::ParameterList& plist);
  ColumnAverageTempEvaluator(const ColumnAverageTempEvaluator& other) = default;
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
  Key temp_key_;
  Key domain_;
  int ncells_depth_;
  double depth_;
private:
  static Utils::RegisteredFactory<Evaluator,ColumnAverageTempEvaluator> reg_;

};
  
} //namespace
} //namespace 

#endif
