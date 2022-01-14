/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the unfrozen fraction model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_SURFACE_KR_EVALUATOR_
#define AMANZI_FLOWRELATIONS_SURFACE_KR_EVALUATOR_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "surface_relperm_model.hh"

namespace Amanzi {
namespace Flow {

class SurfaceRelPermModel;

class SurfaceRelPermEvaluator : public EvaluatorSecondaryMonotypeCV {

 public:
  SurfaceRelPermEvaluator(Teuchos::ParameterList& plist);
  SurfaceRelPermEvaluator(const SurfaceRelPermEvaluator& other);
  Teuchos::RCP<Evaluator> Clone() const;

  Teuchos::RCP<SurfaceRelPermModel> get_Model() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Teuchos::RCP<SurfaceRelPermModel> model_;
  bool is_temp_;
  Key uf_key_;
  Key h_key_;

private:
  static Utils::RegisteredFactory<Evaluator,SurfaceRelPermEvaluator> fac_;


};

} //namespace
} //namespace

#endif

