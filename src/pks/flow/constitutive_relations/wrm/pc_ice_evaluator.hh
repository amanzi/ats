/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  PCIceEvaluator is the interface between state/data and the model, an EOS.

*/

#ifndef AMANZI_RELATIONS_PC_ICE_EVALUATOR_HH_
#define AMANZI_RELATIONS_PC_ICE_EVALUATOR_HH_

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class PCIceWater;

class PCIceEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  // constructor format for all derived classes
  explicit PCIceEvaluator(Teuchos::ParameterList& plist);
  PCIceEvaluator(const PCIceEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;


  Teuchos::RCP<PCIceWater> get_PCIceWater() { return model_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

 protected:
  // the actual model
  Teuchos::RCP<PCIceWater> model_;

  // Keys for fields
  // dependencies
  Key temp_key_;
  Key dens_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, PCIceEvaluator> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
