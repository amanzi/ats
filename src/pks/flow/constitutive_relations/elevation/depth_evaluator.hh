/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#ifndef AMANZI_FLOW_DEPTH_EVALUATOR_HH_
#define AMANZI_FLOW_DEPTH_EVALUATOR_HH_

#include "Factory.hh"
#include "EvaluatorIndependent.hh"

namespace Amanzi {
namespace Flow {

class DepthEvaluator : public EvaluatorIndependentCV {
 public:
  explicit DepthEvaluator(Teuchos::ParameterList& plist);
  DepthEvaluator(const DepthEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from IndependentVariableEvaluator
  virtual void Update_(State& S) override;

 private:
  static Utils::RegisteredFactory<Evaluator, DepthEvaluator> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
