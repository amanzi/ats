/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

A patch evaluator that implements a ponded depth Dirichlet BC

This evaluator is typically used for providing data that are functions of space
and time.  The evaluator consists of a list of region,function pairs, and the
functions are evaluated across that region at each timestep.  If the problem is
time-independent, the `"constant in time`" option results in a performance
boost (as the functions need only be evaluated once).  This leverages the
exaustive functional format capability provided in Amanzi's Functions_ library.

It evaluates into a MultiPatch object

This evaluator is used by providing the option:

`"evaluator type`" == `"flow BC ponded depth`"

.. _flow-bc-ponded-depth-evaluator-spec:
.. admonition:: flow-bc-ponded-depth-evaluator-spec
   KEYS

   - elevation

   INCLUDES

   - independent-variable-patch-function-evaluator-spec

*/

#pragma once

#include "EvaluatorSecondary.hh"
#include "Factory.hh"

namespace Amanzi {

namespace Functions {
class MeshFunction;
}

namespace Flow {
namespace Relations {

class EvaluatorBCPondedDepth
  : public EvaluatorSecondary {

 public:
  explicit EvaluatorBCPondedDepth(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  EvaluatorBCPondedDepth(const EvaluatorBCPondedDepth& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  static const std::string eval_type;
  virtual std::string getType() const override { return eval_type; }

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override {
    // no derivatives make sense here
    return false;
  }

  virtual void EnsureCompatibility(State& S) override;

 protected:
  virtual void Update_(State& S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override {
    AMANZI_ASSERT(false);
  }

 protected:
  Teuchos::RCP<Functions::MeshFunction> func_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorBCPondedDepth> reg_;
};



} // namespace Relations
} // namespace Flow
} // namespace Amanzi
