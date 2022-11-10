#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace ELMKernels {

class CanopyHydrologyEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit CanopyHydrologyEvaluator(Teuchos::ParameterList& plist);
  CanopyHydrologyEvaluator(const CanopyHydrologyEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new CanopyHydrologyEvaluator(*this));
  }

 protected:
  // custom EC used to set subfield names
  virtual void EnsureCompatibility_Structure_(State& S) override;

  // custom EC used because deps have 1 component not 3
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& results) override;

 protected:
  Key domain_;
  Key domain_snow_;

  Key albedo_key_, emissivity_key_;
  Key snow_dens_key_;
  Key unfrozen_fraction_key_;

  double a_water_, a_ice_;
  double e_water_, e_ice_, e_snow_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!

 private:
  static Utils::RegisteredFactory<Evaluator,CanopyHydrologyEvaluator> reg_;
};

} // namespace Relations
} // namespace ELMKernels
