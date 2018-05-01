/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Specifies a value on the surface from the value in the cell just below the
  surface.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_RELATIONS_SURFACE_TOP_CELLS_EVALUATOR_
#define AMANZI_RELATIONS_SURFACE_TOP_CELLS_EVALUATOR_

#include "factory.hh"

#include "EvaluatorSecondary.hh"

namespace Amanzi {
namespace Relations {

class SurfaceTopCellsEvaluator : public EvaluatorSecondary<CompositeVector, CompositeVectorSpace> {

 public:
  explicit
  SurfaceTopCellsEvaluator(Teuchos::ParameterList& plist);

  SurfaceTopCellsEvaluator(const SurfaceTopCellsEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;

  // Required methods from SecondaryVariableFieldEvaluator
  virtual void Evaluate_(const State& S,
                         CompositeVector& result) override;
  
  virtual void EvaluatePartialDerivative_(const State& S, const Key& wrt_key, const Key& wrt_tag,
                                          CompositeVector& result) override {
    ASSERT(0);
  }

  virtual void EnsureCompatibility(State& S) override;

 protected:
  Key dependency_key_, dependency_tag_key_;

 private:
  static Utils::RegisteredFactory<Evaluator,SurfaceTopCellsEvaluator> reg_;

};

} //namespace
} //namespace

#endif
