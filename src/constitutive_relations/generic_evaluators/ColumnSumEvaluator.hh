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
#include "secondary_variable_field_evaluator.hh"

namespace Amanzi {
namespace Relations {

class ColumnSumEvaluator : public SecondaryVariableFieldEvaluator {

public:
  explicit
  ColumnSumEvaluator(Teuchos::ParameterList& plist);
  ColumnSumEvaluator(const ColumnSumEvaluator& other);
  Teuchos::RCP<FieldEvaluator> Clone() const;

protected:
  // Required methods from SecondaryVariableFieldEvaluator
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
               Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);
  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S,
                                      Key request);
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
  static Utils::RegisteredFactory<FieldEvaluator,ColumnSumEvaluator> factory_;

};

} //namespace
} //namespace
