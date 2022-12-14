/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! MeshedElevationEvaluator: evaluates the elevation (z-coordinate) and slope magnitude of a mesh.
/*!
Evaluator type: `"meshed elevation`"

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h z|``

* `"elevation key`" ``[string]`` **elevation** Name the elevation variable. [m]
* `"slope magnitude key`" ``[string]`` **slope_magnitude** Name the elevation variable. [-]
* `"dynamic mesh`" ``[bool]`` **false** Lets the evaluator know that the
  elevation changes in time, and adds the `"deformation`" dependency.
* `"parent domain name`" ``[string]`` **DOMAIN** Domain name of the parent
  mesh, which is the 3D version of this domain.  Attempts to generate an
  intelligent default by stripping "surface" from this domain.

Example:

.. code-block:: xml

  <ParameterList name="elevation">
    <Parameter name="evaluator type" type="string" value="meshed elevation"/>
  </ParameterList>

*/

#ifndef AMANZI_FLOWRELATIONS_MESHED_ELEVATION_EVALUATOR_
#define AMANZI_FLOWRELATIONS_MESHED_ELEVATION_EVALUATOR_

#include "Factory.hh"
#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

namespace Impl {
void
slope_aspect(const AmanziGeometry::Point& normal, double& slope, double& aspect);
}


class MeshedElevationEvaluator : public ElevationEvaluator {
 public:
  explicit MeshedElevationEvaluator(Teuchos::ParameterList& plist);
  MeshedElevationEvaluator(const MeshedElevationEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  virtual void
  EvaluateElevationAndSlope_(const State& S, const std::vector<CompositeVector*>& results) override;

 private:
  static Utils::RegisteredFactory<Evaluator, MeshedElevationEvaluator> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
