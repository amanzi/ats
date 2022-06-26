/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  An elevation evaluator getting values from the volumetric mesh.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_STANDALONE_ELEVATION_EVALUATOR_
#define AMANZI_FLOWRELATIONS_STANDALONE_ELEVATION_EVALUATOR_

#include "Factory.hh"
#include "CompositeVectorFunction.hh"
#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

class StandaloneElevationEvaluator : public ElevationEvaluator {

 public:
  StandaloneElevationEvaluator(Teuchos::ParameterList& elev_plist);
  StandaloneElevationEvaluator(const StandaloneElevationEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void EvaluateElevationAndSlope_(const State& S,
          const std::vector<CompositeVector*>& results) override;

 protected:
  Teuchos::RCP<Functions::CompositeVectorFunction> elevation_function_;
  Teuchos::RCP<Functions::CompositeVectorFunction> slope_function_;
  Teuchos::RCP<Functions::CompositeVectorFunction> aspect_function_;

 private:
  static Utils::RegisteredFactory<Evaluator,StandaloneElevationEvaluator> reg_;
};

} //namespace
} //namespace

#endif
