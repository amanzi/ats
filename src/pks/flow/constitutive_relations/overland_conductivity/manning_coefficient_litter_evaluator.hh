/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The manning coefficient with variable litter evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
Manning's coefficient that varies based on litter thickness and ponded depth.

*/

/*!

(missing documentation!)

*/
#ifndef AMANZI_FLOW_MANNING_COEFFICIENT_LITTER_EVALUATOR_HH_
#define AMANZI_FLOW_MANNING_COEFFICIENT_LITTER_EVALUATOR_HH_

#include "Teuchos_RCP.hpp"

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ManningCoefficientLitterModel;

typedef std::vector<Teuchos::RCP<ManningCoefficientLitterModel>> ManningCoefList;
typedef std::pair<Teuchos::RCP<Functions::MeshPartition>, ManningCoefList> ManningCoefPartition;


class ManningCoefficientLitterEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit ManningCoefficientLitterEvaluator(Teuchos::ParameterList& plist);
  ManningCoefficientLitterEvaluator(const ManningCoefficientLitterEvaluator& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<ManningCoefPartition> get_models() { return models_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  void InitializeFromPlist_();

 protected:
  Key ld_key_;
  Key pd_key_;
  Teuchos::RCP<ManningCoefPartition> models_;

 private:
  static Utils::RegisteredFactory<Evaluator, ManningCoefficientLitterEvaluator> reg_;
};


// Non-member factory
Teuchos::RCP<ManningCoefPartition>
createManningCoefPartition(Teuchos::ParameterList& plist);


} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
