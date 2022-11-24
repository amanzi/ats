/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The IEM Evaluator simply calls the IEM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_ENERGY_RELATIONS_IEM_WATER_VAPOR_EVALUATOR_
#define AMANZI_ENERGY_RELATIONS_IEM_WATER_VAPOR_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "iem_water_vapor.hh"

namespace Amanzi {
namespace Energy {

class IEMWaterVaporEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit IEMWaterVaporEvaluator(Teuchos::ParameterList& plist);
  IEMWaterVaporEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEMWaterVapor>& iem);
  IEMWaterVaporEvaluator(const IEMWaterVaporEvaluator& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<IEMWaterVapor> get_IEM() { return iem_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  void InitializeFromPlist_();

  Key temp_key_;
  Key mol_frac_key_;
  Teuchos::RCP<IEMWaterVapor> iem_;

 private:
  static Utils::RegisteredFactory<Evaluator, IEMWaterVaporEvaluator> factory_;
};

} // namespace Energy
} // namespace Amanzi

#endif
