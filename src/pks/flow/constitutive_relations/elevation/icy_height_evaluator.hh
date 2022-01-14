/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_RELATIONS_ICY_HEIGHT_EVALUATOR_
#define AMANZI_FLOW_RELATIONS_ICY_HEIGHT_EVALUATOR_

#include "height_evaluator.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class IcyHeightModel;

class IcyHeightEvaluator : public HeightEvaluator {

 public:
  // constructor format for all derived classes
  explicit
  IcyHeightEvaluator(Teuchos::ParameterList& plist);
  IcyHeightEvaluator(const IcyHeightEvaluator& other);

  virtual Teuchos::RCP<Evaluator> Clone() const;

  Teuchos::RCP<IcyHeightModel> get_IcyModel() { return icy_model_; }

 protected:
  void InitializeFromPlist_();

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& result);
  virtual void EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
          Key wrt_key, const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key dens_ice_key_;
  Key unfrozen_frac_key_;
  Teuchos::RCP<IcyHeightModel> icy_model_;

 private:
  static Utils::RegisteredFactory<Evaluator,IcyHeightEvaluator> factory_;

};

} //namespace
} //namespace

#endif
