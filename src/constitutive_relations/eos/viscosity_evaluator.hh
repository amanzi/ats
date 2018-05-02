/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSFieldEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_VISC_EVALUATOR_HH_
#define AMANZI_RELATIONS_VISC_EVALUATOR_HH_

#include "viscosity_relation.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Relations {

class ViscosityEvaluator : public EvaluatorSecondary<CompositeVector, CompositeVectorSpace> {

 public:

  // constructor format for all derived classes
  explicit
  ViscosityEvaluator(Teuchos::ParameterList& plist);

  ViscosityEvaluator(const ViscosityEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;


 protected:
  // the actual model
  Teuchos::RCP<ViscosityRelation> visc_;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void Evaluate_(const State& S,
                         CompositeVector& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key, const Key& wrt_tag, CompositeVector& result) override;
  
  // Keys for fields
  // dependencies
  Key temp_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,ViscosityEvaluator> factory_;

};

} // namespace
} // namespace

#endif
