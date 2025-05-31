/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! Computes depth, positive downward relative to the surface, of mesh cells.
/*!

Computes cell depths only.

Note that two algorithms for computing the depth are available:

`"cell centroid`" uses the cell centroid, as computed by the Mesh, directly.

`"mean face centroid`" is the default, as the cell centroid can have problems
in the case of non-planar faces.  Instead, it uses the mean of the above and
below face centroids in place of the face centroid, with the implicit
assumption that dz is uniform (e.g. this is an extruded mesh).

`"evaluator type`" = `"depth`"

.. _evaluator-depth-spec:
.. admonition:: evaluator-depth-spec

   * `"algorithm`" ``[string]`` **mean face centroid** Valid is `"mean face centroid`"
     and `"cell centroid`", see above.

 */


#pragma once

#include "Factory.hh"
#include "EvaluatorIndependent.hh"

namespace Amanzi {
namespace Flow {

class DepthEvaluator : public EvaluatorIndependentCV {
 public:
  explicit DepthEvaluator(Teuchos::ParameterList& plist);
  DepthEvaluator(const DepthEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from IndependentVariableEvaluator
  virtual void Update_(State& S) override;

 protected:
  std::string algorithm_;

 private:
  static Utils::RegisteredFactory<Evaluator, DepthEvaluator> reg_;
};

} // namespace Flow
} // namespace Amanzi
