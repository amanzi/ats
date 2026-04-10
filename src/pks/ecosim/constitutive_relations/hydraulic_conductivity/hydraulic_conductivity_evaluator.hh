/*
  The hydraulic conductivity evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Richards water content evaluator: the standard form as a function of liquid saturation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ECOSIM_HYDRAULIC_CONDUCTIVITY_EVALUATOR_HH_
#define AMANZI_ECOSIM_HYDRAULIC_CONDUCTIVITY_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Ecosim {
namespace Relations {

class HydraulicConductivityModel;

class HydraulicConductivityEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit HydraulicConductivityEvaluator(Teuchos::ParameterList& plist);
  HydraulicConductivityEvaluator(const HydraulicConductivityEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<HydraulicConductivityModel> get_model() { return model_; }

 protected:
   // Required methods from EvaluatorSecondaryMonotypeCV
   virtual void Evaluate_(const State& S,
           const std::vector<CompositeVector*>& result) override;
   virtual void EvaluatePartialDerivative_(const State& S,
           const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

  Key k_key_;
  Key rho_key_;
  Key mu_key_;

  Teuchos::RCP<HydraulicConductivityModel> model_;

 private:
  static Utils::RegisteredFactory<Evaluator,HydraulicConductivityEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace

#endif
