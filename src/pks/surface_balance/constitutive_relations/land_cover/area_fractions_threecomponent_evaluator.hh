/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A subgrid model for determining the area fraction of land, water, and snow within a grid cell.
/*!

Uses a simple linear transition to vary between liquid and bare ground, and
another linear transition to vary between snow-covered and not-snow-covered.

Ordering of the area fractions calculated are: [bare ground, water, snow].

`"evaluator type`" = `"area fractions, three components`"

.. _area_fractions_threecomponent_evaluator-spec:
.. admonition:: area_fractions_threecomponent_evaluator-spec:

   * `"minimum fractional area [-]`" ``[double]`` **1.e-5**
      Mimimum area fraction allowed, less than this is rebalanced as zero.

   DEPENDENCIES:

   - `"snow depth`" **DOMAIN_SNOW-depth**
   - `"ponded depth`" **DOMAIN-ponded_depth**

.. note:

   This evaluator also uses the LandCover_ types.  From that struct, it
   requires the value of the following parameters:

   - `"snow transition height [m]`" ``[double]`` **0.02**
      Minimum thickness for specifying the snow gradient.
   - `"water transition height [m]`" ``[double]`` **0.02**
         Minimum thickness for specifying the water gradient.

.. note:

   This evaluator simplifies the situation by assuming constant density.  This
   make it so that ice and water see the same geometry per unit pressure, which
   isn't quite true thanks to density differences.  However, we hypothesize
   that these differences, on the surface (unlike in the subsurface) really
   don't affect the solution.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class AreaFractionsThreeComponentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit AreaFractionsThreeComponentEvaluator(Teuchos::ParameterList& plist);
  AreaFractionsThreeComponentEvaluator(const AreaFractionsThreeComponentEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new AreaFractionsThreeComponentEvaluator(*this));
  }

 protected:
  // custom EC used to set subfield names
  virtual void EnsureCompatibility_Structure_(State& S) override;

  // custom EC used because deps have 1 component not 3
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override
  {
    result[0]->PutScalar(0.);
    //Errors::Message msg("NotImplemented: AreaFractionsThreeComponentEvaluator currently does not provide derivatives.");
    //Exceptions::amanzi_throw(msg);
  }

 protected:
  Key domain_, domain_snow_;
  Key snow_depth_key_;
  Key ponded_depth_key_;
  double min_area_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, AreaFractionsThreeComponentEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
