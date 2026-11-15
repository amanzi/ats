#ifndef ATS_ECOSIM_LAYER0_EVALUATOR_HH_
#define ATS_ECOSIM_LAYER0_EVALUATOR_HH_
#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"
namespace Amanzi { namespace EcoSIM {
// Add retained native litter water and dry litter heat to the existing mobile surface stores.
class EcoLayer0Evaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit EcoLayer0Evaluator(Teuchos::ParameterList& plist);
  Teuchos::RCP<Evaluator> Clone() const override { return Teuchos::rcp(new EcoLayer0Evaluator(*this)); }
 protected:
  void EnsureCompatibility_ToDeps_(State& S) override;
  void Evaluate_(const State& S,const std::vector<CompositeVector*>& result) override;
  void EvaluatePartialDerivative_(const State& S,const Key& wrt,const Tag& tag,
                                  const std::vector<CompositeVector*>& result) override;
 private:
  void Calculate_(const State& S,const Key& wrt,const std::vector<CompositeVector*>& result) const;
  bool energy_;
  Key base_,material_,pressure_,area_,temperature_,fraction_,liquid_,ice_;
  static Utils::RegisteredFactory<Evaluator,EcoLayer0Evaluator> reg_;
};
} }
#endif
