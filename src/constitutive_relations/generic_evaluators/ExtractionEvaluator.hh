/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

Extracts a field on one mesh from a field on a superset of that mesh using
parent entities.

*/
#pragma once

#include "Factory.hh"

#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Relations {

class ExtractionEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit ExtractionEvaluator(Teuchos::ParameterList& plist);
  ExtractionEvaluator(const ExtractionEvaluator& other) = default;
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
    AMANZI_ASSERT(false);
  }

  Key dependency_key_;
  Key parent_domain_;

 private:
  static Utils::RegisteredFactory<Evaluator, ExtractionEvaluator> reg_;
};

} // namespace Relations
} // namespace Amanzi
