/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Evaluator for determining height( rho, head )

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
  explicit IcyHeightEvaluator(Teuchos::ParameterList& plist);
  IcyHeightEvaluator(const IcyHeightEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<IcyHeightModel> get_IcyModel() { return icy_model_; }

 protected:
  void InitializeFromPlist_();

  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  Key dens_ice_key_;
  Key unfrozen_frac_key_;
  Teuchos::RCP<IcyHeightModel> icy_model_;

 private:
  static Utils::RegisteredFactory<Evaluator, IcyHeightEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
