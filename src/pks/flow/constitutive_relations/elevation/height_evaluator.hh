/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_HEIGHT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_HEIGHT_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class HeightModel;

class HeightEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit HeightEvaluator(Teuchos::ParameterList& plist);
  HeightEvaluator(const HeightEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<HeightModel> get_Model() { return model_; }

  void set_bar(bool bar) { bar_ = bar; }

 protected:
  // Needs a special EnsureCompatibility to get around trying to find face
  // values and derivatives of face values.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key dens_key_;
  Key pres_key_;
  Key gravity_key_;
  Key patm_key_;
  bool bar_;

  Teuchos::RCP<HeightModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator, HeightEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
