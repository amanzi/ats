/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

*/

#ifndef AMANZI_RELATIONS_OVERLAND_SOURCE_FROM_SUBSURFACE_FLUX_EVALUATOR_HH_
#define AMANZI_RELATIONS_OVERLAND_SOURCE_FROM_SUBSURFACE_FLUX_EVALUATOR_HH_

#include "Evaluator_Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class OverlandSourceFromSubsurfaceFluxEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit OverlandSourceFromSubsurfaceFluxEvaluator(Teuchos::ParameterList& plist);
  OverlandSourceFromSubsurfaceFluxEvaluator(
    const OverlandSourceFromSubsurfaceFluxEvaluator& other) = default;
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

  void IdentifyFaceAndDirection_(const State& S);

  typedef std::pair<int, double> FaceDir;
  Teuchos::RCP<std::vector<FaceDir>> face_and_dirs_;

  Key flux_key_;
  Key dens_key_;
  bool volume_basis_;

  Key domain_surf_;
  Key domain_sub_;

 private:
  static Utils::RegisteredFactory<Evaluator, OverlandSourceFromSubsurfaceFluxEvaluator> fac_;
};

} // namespace Relations
} // namespace Amanzi

#endif
