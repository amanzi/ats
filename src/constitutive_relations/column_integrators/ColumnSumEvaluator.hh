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
     variable's domain (typically "surface" or "surface_column:\*") and is
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
  IntegratorColumnSum(Teuchos::ParameterList& plist,
                      std::vector<const Epetra_MultiVector*>& deps,
                      const AmanziMesh::Mesh* mesh);
  int scan(AmanziMesh::Entity_ID col, AmanziMesh::Entity_ID c, AmanziGeometry::Point& p);
  double coefficient(AmanziMesh::Entity_ID col);

 private:
  bool volume_average_;
  bool volume_factor_;
  bool divide_by_density_;
  double coef_;
  const Epetra_MultiVector* integrand_;
  const Epetra_MultiVector* cv_;
  const Epetra_MultiVector* surf_cv_;
  const Epetra_MultiVector* dens_;
};

} // namespace Impl

using ColumnSumEvaluator =
  EvaluatorColumnIntegrator<Impl::ParserColumnSum, Impl::IntegratorColumnSum>;


} // namespace Relations
} // namespace Amanzi
