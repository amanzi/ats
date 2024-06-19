/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
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

#include <functional>
#include "Factory.hh"
#include "BlockVector_decl.hh"
#include "EvaluatorColumnIntegrator.hh"

namespace Amanzi {
namespace Relations {

namespace Impl {

struct ParserColumnSum {
  ParserColumnSum(Teuchos::ParameterList& plist, const KeyTag& key_tag);
  KeyTagSet dependencies;
};


class IntegratorColumnSum {
 public:
  using cView_type = BlockVector<double>::cView_type;

  IntegratorColumnSum(Teuchos::ParameterList& plist,
                      std::vector<cView_type>& deps,
                      const AmanziMesh::Mesh& mesh);

  KOKKOS_INLINE_FUNCTION
  int scan(const AmanziMesh::Entity_ID col,
           const AmanziMesh::Entity_ID c,
           AmanziGeometry::Point& p) const;

  KOKKOS_INLINE_FUNCTION
  double coefficient(const AmanziMesh::Entity_ID col) const;

 private:
  bool volume_average_;
  bool volume_factor_;
  bool divide_by_density_;
  double coef_;
  cView_type integrand_;
  cView_type cv_;
  cView_type surf_cv_;
  cView_type dens_;
};


KOKKOS_INLINE_FUNCTION
int
IntegratorColumnSum::scan(const AmanziMesh::Entity_ID col,
                          const AmanziMesh::Entity_ID c,
                          AmanziGeometry::Point& p) const
{
  double contrib = integrand_(c, 0);
  if (volume_average_ || volume_factor_) { contrib *= cv_(c, 0); }
  if (divide_by_density_) { contrib /= dens_(c, 0); }
  p[0] += contrib;

  if (volume_average_) p[1] += cv_(c, 0);
  return false;
}


KOKKOS_INLINE_FUNCTION
double
IntegratorColumnSum::coefficient(const AmanziMesh::Entity_ID col) const
{
  if (volume_factor_) {
    return coef_ / surf_cv_(col, 0);
  } else {
    return coef_;
  }
}


} // namespace Impl

using ColumnSumEvaluator =
  EvaluatorColumnIntegrator<Impl::ParserColumnSum, Impl::IntegratorColumnSum>;

} // namespace Relations
} // namespace Amanzi
