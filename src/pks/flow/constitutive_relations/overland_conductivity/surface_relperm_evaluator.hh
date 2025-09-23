/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the unfrozen fraction model.

*/

#ifndef AMANZI_FLOWRELATIONS_SURFACE_KR_EVALUATOR_
#define AMANZI_FLOWRELATIONS_SURFACE_KR_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "surface_relperm_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class SurfaceRelPermModel;

class SurfaceRelPermEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  SurfaceRelPermEvaluator(Teuchos::ParameterList& plist);
  SurfaceRelPermEvaluator(const SurfaceRelPermEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<SurfaceRelPermModel> get_Model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Teuchos::RCP<SurfaceRelPermModel> model_;
  bool is_temp_;
  Key uf_key_;
  Key h_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, SurfaceRelPermEvaluator> fac_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
