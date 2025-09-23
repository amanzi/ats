/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky  (dasvyat@lanl.gov)
*/

/*
  An evaluator for converting the darcy flux to volumetric flux

*/
#ifndef AMANZI_RELATIONS_VOL_DARCY_FLUX_HH_
#define AMANZI_RELATIONS_VOL_DARCY_FLUX_HH_

#include "Evaluator_Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ATS_Physics {
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
} // namespace ATS_Physics
} // namespace Amanzi

#endif
