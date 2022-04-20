/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A subgrid model for determining the area fraction of snow vs not snow within a grid cell.
/*!

Uses a simple linear transition to vary between liquid and bare ground, and
another linear transition to vary between snow-covered and not-snow-covered.

Ordering of the area fractions calculated are: [bare ground/water, snow].

.. _area-fractions-twocomponent-evaluator-spec:
.. admonition:: area-fractions-twocomponent-evaluator-spec:

   * `"minimum fractional area [-]`" ``[double]`` **1.e-5**
         Mimimum area fraction allowed, less than this is rebalanced as zero.

   KEYS:
   DEPENDENCIES:

   - `"snow depth`" ``[string]``

.. note:

   This evaluator also uses the LandCover_ types.  From that struct, it
   requires the value of the following parameters:

   - `"snow transition height [m]`" ``[double]`` **0.02**
      Minimum thickness for specifying the snow gradient.

.. note:

   This evaluator simplifies the situation by assuming constant density.  This
   make it so that ice and water see the same geometry per unit pressure, which
   isn't quite true thanks to density differences.  However, we hypothesize
   that these differences, on the surface (unlike in the subsurface) really
   don't affect the solution.

*/

#pragma once

#include "Factory.hh"
#include "secondary_variable_field_evaluator.hh"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class AreaFractionsTwoComponentEvaluator : public SecondaryVariableFieldEvaluator {

 public:
  explicit
  AreaFractionsTwoComponentEvaluator(Teuchos::ParameterList& plist);
  AreaFractionsTwoComponentEvaluator(const AreaFractionsTwoComponentEvaluator& other) = default;

  virtual Teuchos::RCP<FieldEvaluator> Clone() const override {
    return Teuchos::rcp(new AreaFractionsTwoComponentEvaluator(*this));
  }

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S) override;

 protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result) override;
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) override;

 protected:
  Key domain_, domain_snow_;
  Key snow_depth_key_;
  double min_area_;

  // this is horrid, because this cannot yet live in state
  // bring on new state!
  LandCoverMap land_cover_;

 private:
  static Utils::RegisteredFactory<FieldEvaluator,AreaFractionsTwoComponentEvaluator> reg_;

};

} //namespace
} //namespace
} //namespace
