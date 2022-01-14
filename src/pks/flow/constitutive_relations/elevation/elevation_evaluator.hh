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
  explicit
  ElevationEvaluator(Teuchos::ParameterList& plist);

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results);

  virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results) = 0;

  virtual bool HasFieldChanged(const Teuchos::Ptr<State>& S, Key request);

  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

protected:
  bool updated_once_;
  bool dynamic_mesh_;
  Key deformation_key_;

};

} //namespace
} //namespace

#endif
