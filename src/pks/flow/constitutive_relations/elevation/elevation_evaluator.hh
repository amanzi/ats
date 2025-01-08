/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The elevation evaluator gets the surface elevation and slope.

  This is not a normal EvaluatorSecondaryMonotypeCV, as it has no
  dependencies, which means we have to force it to update (dependencies
  will never have changed) in HasFieldChanged.  This is done this
  way so that when the mesh changes, this can be updated appropriately.

*/

#ifndef AMANZI_FLOWRELATIONS_ELEVATION_EVALUATOR_
#define AMANZI_FLOWRELATIONS_ELEVATION_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {

class ElevationEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit ElevationEvaluator(Teuchos::ParameterList& plist);
  ElevationEvaluator(const ElevationEvaluator& other) = default;

  virtual bool Update(State& S, const Key& request) override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  virtual void
  EvaluateElevationAndSlope_(const State& S, const std::vector<CompositeVector*>& results) = 0;

  //
  // This is required to make sure that elevation, slope, and aspect share a
  // common structure.  Often, aspect is not used and so it can otherwise be an
  // empty vector with no structure, which causes seg faults.
  //
  virtual void EnsureCompatibility_Structure_(State& S) override
  {
    EnsureCompatibility_StructureSame_(S);
  }


 protected:
  bool updated_once_;
  bool dynamic_mesh_;
  Key deformation_key_;
};

} // namespace Flow
} // namespace Amanzi

#endif
