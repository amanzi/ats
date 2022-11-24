/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_VISC_EVALUATOR_HH_
#define AMANZI_RELATIONS_VISC_EVALUATOR_HH_

#include "viscosity_relation.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class ViscosityEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit ViscosityEvaluator(Teuchos::ParameterList& plist);

  ViscosityEvaluator(const ViscosityEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  // the actual model
  Teuchos::RCP<ViscosityRelation> visc_;

  // Keys for fields
  // dependencies
  Key temp_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, ViscosityEvaluator> factory_;
};

} // namespace Relations
} // namespace Amanzi

#endif
