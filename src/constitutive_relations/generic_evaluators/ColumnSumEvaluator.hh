/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/
//! Sums a subsurface field vertically only a surface field.

/*!

Simple vertical sum of all cells below each surface cell.  Note that their are
options for including volume factors (multiply by subsurface cell volume, sum,
divide by surface cell area) and density (useful for the most common use case
of summing fluxes onto the surface and converting to m/s instead of mol/m^2/s).


.. _column-sum-evaluator-spec:
.. admonition:: column-sum-evaluator-spec

    * `"include volume factor`" ``[bool]`` **true** In summing, multiply the
      summand subsurface cell volume, then divide the sum by the surface cell
      area.  Useful for converting volumetric fluxes to total fluxes.

    * `"divide by density`" ``[bool]`` **true** Divide the summand by density.
      Useful for converting molar fluxes to volumetric fluxes
      (e.g. transpiration).

    * `"column domain name`" ``[string]`` **domain** The domain of the
      subsurface mesh.  Note this defaults to a sane thing based on the
      variable's domain (typically "surface" or "surface_column:*") and is
      rarely set by the user.

    KEYS:
    - `"summed`" The summand, defaults to the root suffix of the calculated variable.
    - `"cell volume`" Defaults to domain's cell volume.
    - `"surface cell volume`" Defaults to surface domain's cell volume.
    - `"molar density`" Defaults to domain's molar_density_liquid.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class ColumnSumEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  explicit
  ColumnSumEvaluator(Teuchos::ParameterList& plist);
  ColumnSumEvaluator(const ColumnSumEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  virtual void EnsureCompatibility(State& S) override;
  virtual bool Update(State& S, const Key& request) override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
                              const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
               const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result) override;

 protected:
  double coef_;

  Key dep_key_;
  Key cv_key_;
  Key molar_dens_key_;
  Key surf_cv_key_;

  Key domain_;
  Key surf_domain_;

  bool updated_once_;
private:
  static Utils::RegisteredFactory<Evaluator,ColumnSumEvaluator> factory_;

};

} //namespace
} //namespace
