/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A subgrid model for determining the area fraction of land, open water, and snow within a grid cell.
/*!

Uses the subgrid equation from Jan et al WRR 2018 for volumetric or effective
ponded depth to determine the area of water, then heuristically places snow on
top of that surface.

`"evaluator type`" = `"area fractions, three components with microtopography`"

.. _evaluator-area-fractions-three-components-with-microtopography-spec:
.. admonition:: evaluator-area-fractions-three-components-with-microtopography-spec

   * `"snow transitional height [m]`" ``[double]`` **0.02**
     Minimum thickness for specifying the snow gradient.
   * `"minimum fractional area [-]`" ``[double]`` **1.e-5**
     Mimimum area fraction allowed, less than this is rebalanced as zero.
   * `"snow domain name`" ``[string]`` **DOMAIN_SNOW** A default is guessed at
     by replacing `"surface`" with `"snow`" in the this's domain.

   KEYS:

   - `"microtopographic relief`" **DOMAIN-microtopographic_relief**
     The name of del_max, the max microtopography value.
   - `"excluded volume`" **DOMAIN-excluded_volume**
     The name of del_excluded, the integral of the microtopography.
   - `"ponded depth`" **DOMAIN-pressure**
     The name of the surface water ponded depth.
   - `"snow depth`" **DOMAIN_SNOW-depth**
     The name of the snow depth.
   - `"volumetric snow depth`" **DOMAIN_SNOW-volumetric_depth**
     The name of the snow depth.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class AreaFractionsThreeComponentMicrotopographyEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit AreaFractionsThreeComponentMicrotopographyEvaluator(Teuchos::ParameterList& plist);
  AreaFractionsThreeComponentMicrotopographyEvaluator(
    const AreaFractionsThreeComponentMicrotopographyEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new AreaFractionsThreeComponentMicrotopographyEvaluator(*this));
  }

 protected:
  // custom EC used to set subfield names
  void EnsureCompatibility_Structure_(State& S) override;

  // custom EC used because deps have 1 component not 3
  void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override
  {
    Exceptions::amanzi_throw("NotImplemented: AreaFractionsThreeComponentMicrotopographyEvaluator "
                             "currently does not provide derivatives.");
  }

 protected:
  Key domain_, domain_snow_;
  Key ponded_depth_key_, snow_depth_key_, vol_snow_depth_key_;
  Key delta_max_key_, delta_ex_key_;
  double rho_liq_;
  double snow_subgrid_transition_;
  double min_area_;

 private:
  static Utils::RegisteredFactory<Evaluator, AreaFractionsThreeComponentMicrotopographyEvaluator>
    reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
