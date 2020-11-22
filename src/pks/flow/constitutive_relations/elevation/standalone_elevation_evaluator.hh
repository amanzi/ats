/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  An elevation evaluator getting values from the volumetric mesh.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_STANDALONE_ELEVATION_EVALUATOR_
#define AMANZI_FLOWRELATIONS_STANDALONE_ELEVATION_EVALUATOR_

#include "CompositeVectorFunction.hh"
#include "elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

class StandaloneElevationEvaluator : public ElevationEvaluator {

 public:
  StandaloneElevationEvaluator(Teuchos::ParameterList& elev_plist);
  StandaloneElevationEvaluator(const StandaloneElevationEvaluator& other);

  Teuchos::RCP<FieldEvaluator> Clone() const;

  virtual void EvaluateElevationAndSlope_(const Teuchos::Ptr<State>& S,
          const std::vector<Teuchos::Ptr<CompositeVector> >& results);

 protected:
  Teuchos::RCP<Functions::CompositeVectorFunction> elevation_function_;
  Teuchos::RCP<Functions::CompositeVectorFunction> slope_function_;
  Teuchos::RCP<Functions::CompositeVectorFunction> aspect_function_;

};

} //namespace
} //namespace

#endif
