/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Provides a depth-based profile of root density.
/*!

Sets the (discrete) root fraction as a function of depth.  The rooting density
is given by:

.. math:
   \rho_root =  \frac{1}{2} ( \alpha \; exp(-\alpha z) + \beta \; exp(-\beta z) )

This function is such that the integral over depth = [0,inf) is 1.  Then,
computing this over the vertical corridor is done by integrating this function
between the depth of the face above and below for each grid cell, with the
bottom-most grid cell integrating to infinity.

Note that the two parameters, :math:`\alpha` and :math:`\beta` are provided in
the Land Cover class as `"rooting profile alpha`" and `"rooting profile beta`"
respectively.

.. _rooting-depth-fraction-evaluator-spec:
.. admonition:: rooting-depth-fraction-evaluator-spec

   * `"surface domain name`" ``[string]`` **SURFACE_DOMAIN** Sane default provided for most domain names.

   KEYS:

   - `"cell volume`" **DOMAIN-cell_volume**
   - `"surface area`" **SURFACE_DOMAIN-cell_volume**

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class RootingDepthFractionEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit RootingDepthFractionEvaluator(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  RootingDepthFractionEvaluator(const RootingDepthFractionEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return false;
  }

  // need a custom EnsureCompatibility as some vectors cross meshes.
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  void InitializeFromPlist_();

  std::string getType() const override { return eval_type; }

 protected:
  static const std::string eval_type;

  Key cv_key_;
  Key surf_cv_key_;

  Key domain_surf_;
  Key domain_sub_;

  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, RootingDepthFractionEvaluator> reg_;
};


namespace Impl {

KOKKOS_INLINE_FUNCTION
double computeIntegralRootFunc(double z, double alpha, double beta) {
  return -0.5 * (Kokkos::exp(-alpha * z) + Kokkos::exp(-beta * z));
}

}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
