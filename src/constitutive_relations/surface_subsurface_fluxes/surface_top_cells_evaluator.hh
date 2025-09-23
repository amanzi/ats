/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Specifies a value on the surface from the value in the cell just below the
  surface.

*/

#ifndef AMANZI_RELATIONS_SURFACE_TOP_CELLS_EVALUATOR_
#define AMANZI_RELATIONS_SURFACE_TOP_CELLS_EVALUATOR_

#include "Factory.hh"

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

class SurfaceTopCellsEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SurfaceTopCellsEvaluator(Teuchos::ParameterList& plist);
  SurfaceTopCellsEvaluator(const SurfaceTopCellsEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override;

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override
  {
    AMANZI_ASSERT(0);
  }

 protected:
  Key dependency_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, SurfaceTopCellsEvaluator> reg_;
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi

#endif
