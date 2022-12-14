/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  An evaluator for converting the darcy flux to volumetric flux

  Authors: Daniil Svyatsky  (dasvyat@lanl.gov)
*/
#ifndef AMANZI_RELATIONS_VOL_DARCY_FLUX_HH_
#define AMANZI_RELATIONS_VOL_DARCY_FLUX_HH_

#include "Evaluator_Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class Volumetric_FluxEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit Volumetric_FluxEvaluator(Teuchos::ParameterList& plist);
  Volumetric_FluxEvaluator(const Volumetric_FluxEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // custom ensure compatibility as all data is not just on the same components
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key flux_key_;
  Key dens_key_;
  Key mesh_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, Volumetric_FluxEvaluator> fac_;
};


} // namespace Relations
} // namespace Amanzi

#endif
