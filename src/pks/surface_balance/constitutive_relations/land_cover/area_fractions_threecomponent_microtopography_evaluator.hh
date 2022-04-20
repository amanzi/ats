/*
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

.. _area-fractions-threecomponent-microtopography-evaluator-spec:
.. admonition:: area-fractions-threecomponent-microtopography-evaluator-spec

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


NOTE: this evaluator simplifies the situation by assuming constant density.
This make it so that ice and water see the same geometry per unit pressure,
which isn't quite true thanks to density differences.  However, we hypothesize
that these differences, on the surface (unlike in the subsurface) really don't
matter much. --etc

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class AreaFractionsThreeComponentMicrotopographyEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  AreaFractionsThreeComponentMicrotopographyEvaluator(Teuchos::ParameterList& plist);
  AreaFractionsThreeComponentMicrotopographyEvaluator(const AreaFractionsThreeComponentMicrotopographyEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const {
    return Teuchos::rcp(new AreaFractionsThreeComponentMicrotopographyEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {
    Exceptions::amanzi_throw("NotImplemented: AreaFractionsThreeComponentMicrotopographyEvaluator currently does not provide derivatives.");
  }

 protected:

  Key domain_, domain_snow_;
  Key ponded_depth_key_, snow_depth_key_, vol_snow_depth_key_;
  Key delta_max_key_, delta_ex_key_;
  double rho_liq_;
  double snow_subgrid_transition_;
  double min_area_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,AreaFractionsThreeComponentMicrotopographyEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace
