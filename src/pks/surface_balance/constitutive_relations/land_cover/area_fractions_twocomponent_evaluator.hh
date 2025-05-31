/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A subgrid model for determining the area fraction of snow vs not snow within a grid cell.
/*!

Uses a simple linear transition to vary between liquid and bare ground, and
another linear transition to vary between snow-covered and not-snow-covered.

Ordering of the area fractions calculated are: [bare ground/water, snow].

This evaluator simplifies the situation by assuming constant water density.
This make it so that ice and water see the same geometry per unit pressure,
which isn't quite true thanks to density differences.  However, we hypothesize
that these differences, on the surface (unlike in the subsurface) really don't
affect the solution.


`"evaluator type`" = `"area fractions, two components`"

.. _evaluator-area-fractions-two-components-spec:
.. admonition:: evaluator-area-fractions-two-components-spec:

   * `"minimum fractional area [-]`" ``[double]`` **1.e-5** Minimum area
     fraction allowed, less than this is rebalanced as zero.

   DEPENDENCIES:

   - `"snow depth`" ``[string]``

.. note::

   This evaluator also uses the :ref:`Land Cover` types.  From that struct, it
   requires the value of the following parameters:

   - `"snow transition height [m]`"

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class AreaFractionsTwoComponentEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit AreaFractionsTwoComponentEvaluator(Teuchos::ParameterList& plist);
  AreaFractionsTwoComponentEvaluator(const AreaFractionsTwoComponentEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new AreaFractionsTwoComponentEvaluator(*this));
  }

 protected:
  // custom EC used to set subfield names
  virtual void EnsureCompatibility_Structure_(State& S) override;

  // custom EC used because deps have 1 component not 2
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key domain_, domain_snow_;
  Key snow_depth_key_;
  double min_area_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<Evaluator, AreaFractionsTwoComponentEvaluator> reg_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
