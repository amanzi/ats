/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Provides a depth-based profile of root density.
/*!

Sets the root fraction as a function of depth,

.. math:
   F_root =  ( \alpha \; exp(-\alpha z) + \beta \; exp(-\beta z) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.

Note that all three parameters, a, b, and the cutoff, are provided in the
LandCover type.

.. _rooting-depth-fraction-evaluator-spec:
.. admonition:: rooting-depth-fraction-evaluator-spec

   * `"surface domain name`" ``[string]`` **SURFACE_DOMAIN** Sane default provided for most domain names.

   KEYS:

   - `"depth`" **DOMAIN-depth**
   - `"cell volume`" **DOMAIN-cell_volume**
   - `"surface cell volume`" **SURFACE_DOMAIN-cell_volume**

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class RootingDepthFractionModel;

class RootingDepthFractionEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit RootingDepthFractionEvaluator(Teuchos::ParameterList& plist);
  RootingDepthFractionEvaluator(const RootingDepthFractionEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  void InitializeFromPlist_();

 protected:
  Key z_key_;
  Key cv_key_;
  Key surf_cv_key_;

  Key domain_surf_;
  Key domain_sub_;

  LandCoverMap land_cover_;
  std::map<std::string, Teuchos::RCP<RootingDepthFractionModel>> models_;

 private:
  static Utils::RegisteredFactory<Evaluator, RootingDepthFractionEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
