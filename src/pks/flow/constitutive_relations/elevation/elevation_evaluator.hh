/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The elevation evaluator gets the surface elevation and slope.

  This is not a normal EvaluatorSecondaryMonotypeCV, as it has no
  dependencies, which means we have to force it to update (dependencies
  will never have changed) in HasFieldChanged.  This is done this
  way so that when the mesh changes, this can be updated appropriately.

  Authors: Ethan Coon (ecoon@lanl.gov)
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

  // I don't believe this should need to be implemented. --ETC
  //  virtual void EnsureCompatibility(State& S) override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S,
          const std::vector<CompositeVector*>& results) override;
  virtual void EvaluatePartialDerivative_(const State& S,
          const Key& wrt_key, const Tag& wrt_tag,
          const std::vector<CompositeVector*>& results) override;

  virtual void EvaluateElevationAndSlope_(const State& S,
          const std::vector<CompositeVector*>& results) = 0;

 protected:
  bool updated_once_;
  bool dynamic_mesh_;
  Key deformation_key_;

};

} //namespace
} //namespace

#endif
