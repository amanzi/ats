/*
  Interception/throughfall rate.

  A simple model based on relaxation from current water content to a saturated water content.


  Model:

          |
          | source
          V
         /   \
      I /     \  T
       V       |
  --Theta--    |
       |       |
       | D     | 
       V       V
  ----------------------

  This is the model for splitting source into throughfall and interception.  

  
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_INTERCEPTION_EVALUATOR_HH_
#define AMANZI_RELATIONS_INTERCEPTION_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class InterceptionEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:

  // constructor format for all derived classes
  explicit
  InterceptionEvaluator(Teuchos::ParameterList& plist);

  InterceptionEvaluator(const InterceptionEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key ai_key_;
  Key source_key_;

  double n_liq_;
  double throughfall_coef_;
  bool source_in_meters_;

 private:
  static Amanzi::Utils::RegisteredFactory<Evaluator,InterceptionEvaluator> factory_;

};

} // namespace
} // namespace
} // namespace

#endif
