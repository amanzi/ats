/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! SubgridMobileDepthEvaluator: calculates mobile depth including a depression storage term.
/*!

Effectively, this is

.. math:
   \delta_{mobile} = max(0, \delta - \delta_{depression})

from Jan et al WRR 2018.

It is a bit trickier than this because depression depth is only provided on
cells, while mobile depth (and ponded depth) are on both cells and boundary
faces.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace FlowRelations {

class SubgridMobileDepthEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SubgridMobileDepthEvaluator(Teuchos::ParameterList& plist);
  SubgridMobileDepthEvaluator(const SubgridMobileDepthEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 private:
  Key depth_key_;
  Key depr_depth_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, SubgridMobileDepthEvaluator> factory_;
};

} // namespace FlowRelations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
