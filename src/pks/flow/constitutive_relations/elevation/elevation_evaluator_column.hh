/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! ElevationEvaluatorColumn: evaluates the elevation (z-coordinate) and slope magnitude of a mesh.

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

/*!
Evaluator type: `"elevation column`"

Evaluates the z-coordinate and the magnitude of the slope :math:``|\nambla_h z|``

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
    <Parameter name="field evaluator type" type="string" value="column elevation"/>
  </ParameterList>

*/

#ifndef AMANZI_FLOWRELATIONS_ELEVATION_EVALUATOR_COLUMN_
#define AMANZI_FLOWRELATIONS_ELEVATION_EVALUATOR_COLUMN_

#include "Factory.hh"
#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

class ElevationEvaluatorColumn : public ElevationEvaluator {

 public:
  explicit
  ElevationEvaluatorColumn(Teuchos::ParameterList& plist);
  ElevationEvaluatorColumn(const ElevationEvaluatorColumn& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void EvaluateElevationAndSlope_(const State& S,
          const std::vector<CompositeVector*>& results) override;

  // Custom EnsureCompatibility fills dependencies based on domain set.
  virtual void EnsureCompatibility(State& S) override;

 private:
  static Utils::RegisteredFactory<Evaluator,ElevationEvaluatorColumn> reg_;

  Key base_poro_suffix_;
  Key surface_domain_;
  Key dset_name_;
};

} //namespace
} //namespace

#endif
