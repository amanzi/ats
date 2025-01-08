/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

//! ElevationEvaluatorColumn: evaluates the elevation (z-coordinate) and slope magnitude of a mesh.
/*!

Evaluates elevation, slope, and aspect of the "surface_star" mesh of the Arctic
Intermediate Scale Model (ISM).

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h
z|`` Note that this is evaluator is different from both the standard elevation
evaluator, which would take elevation/slope/aspect from the parent 3D mesh, and
from the standard "single cell" column mesh elevation, which can do the same.
Instead, it is a mix of:

- elevation is aggregated from the column meshes' top faces

- slope is calculated using a cell-neighbor differencing scheme.  Unlike 3D
  meshes, in the ISM, we deform columns of cells and so do not have faces.
  Therefore, we use all neighboring cells to compute an average normal for the
  cell, and use that to compute slope and aspect.

- aspect is set to 0.  It could easily be calculated using the same normal as
  the slope algorithm, but is not done currently.

`"evaluator type`" = `"elevation column`"

.. _column_elevation_evaluator-spec:
.. admonition:: column_elevation_evaluator-spec

   * `"elevation key`" ``[string]`` **elevation** Name the elevation variable. [m]
   * `"slope magnitude key`" ``[string]`` **slope_magnitude** Name the elevation
     variable. [-]
   * `"dynamic mesh`" ``[bool]`` **false** Lets the evaluator know that the
     elevation changes in time, and adds the `"deformation`" and
     `"base_porosity`" dependencies.
   * `"parent domain name`" ``[string]`` **DOMAIN** Domain name of the parent
     mesh, which is the 3D version of this domain.  In the columnar meshes the
     surface elevation and slope are assigned based on the columns and not the
     base 3D domain.

Example:

.. code-block:: xml

  <ParameterList name="surface_star-elevation">
    <Parameter name="evaluator type" type="string" value="column elevation"/>
  </ParameterList>

*/

#ifndef AMANZI_FLOWRELATIONS_ELEVATION_EVALUATOR_COLUMN_
#define AMANZI_FLOWRELATIONS_ELEVATION_EVALUATOR_COLUMN_

#include "Factory.hh"
#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

class ColumnElevationEvaluator : public ElevationEvaluator {
 public:
  explicit ColumnElevationEvaluator(Teuchos::ParameterList& plist);
  ColumnElevationEvaluator(const ColumnElevationEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void
  EvaluateElevationAndSlope_(const State& S, const std::vector<CompositeVector*>& results) override;

  // Custom EnsureCompatibility fills dependencies based on domain set.
  virtual void EnsureEvaluators(State& S) override;

  // do not pass my structure to my dependencies
  virtual void EnsureCompatibility_ToDeps_(State& S) override {}


 private:
  static Utils::RegisteredFactory<Evaluator, ColumnElevationEvaluator> reg_;

  Key base_poro_suffix_;
  Key surface_domain_;
  Key dset_name_;
};

} // namespace Flow
} // namespace Amanzi

#endif
