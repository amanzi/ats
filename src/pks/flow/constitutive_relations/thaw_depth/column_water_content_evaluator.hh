/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The column water content evaluator gets the subsurface water content.
  This computes the water (frozen + unfrozen) content in the column and puts on the surface cell.
  This is EvaluatorSecondaryMonotypeCV and depends on the subsurface water content, 

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_COLWATERCONTENT_EVALUATOR_
#define AMANZI_FLOWRELATIONS_COLWATERCONTENT_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ColumnWaterContentEvaluator : public EvaluatorSecondaryMonotypeCV {

public:
  explicit
  ColumnWaterContentEvaluator(Teuchos::ParameterList& plist);
  ColumnWaterContentEvaluator(const ColumnWaterContentEvaluator& other);
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
  Key  wc_key_;
  Key domain_;
private:
  static Utils::RegisteredFactory<Evaluator,ColumnWaterContentEvaluator> reg_;

};
  
} //namespace
} //namespace 

#endif
