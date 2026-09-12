/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!
A patch evaluator that implements a seepage face head BC for overland flow.

On each face, given a boundary head h_f [m] and the adjacent cell's ponded
depth h_c and elevations z_f, z_c: if h_f + z_f >= h_c + z_c the face is a
no-flow (Neumann 0) boundary, otherwise it is a Dirichlet boundary with water
level h_f + z_f.  This matches the CPU "seepage face head" condition.

`"evaluator type`" == `"flow BC seepage face head`"

.. _flow-bc-seepage-face-head-evaluator-spec:
.. admonition:: flow-bc-seepage-face-head-evaluator-spec

   KEYS
   - elevation
   - ponded depth

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

class EvaluatorBCSeepageFaceHead : public EvaluatorSecondary {
 public:
  explicit EvaluatorBCSeepageFaceHead(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  EvaluatorBCSeepageFaceHead(const EvaluatorBCSeepageFaceHead& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  static const std::string eval_type;
  virtual std::string getType() const override { return eval_type; }

  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    return false;
  }

  virtual void EnsureCompatibility(State& S) override;

  // not protected because it calls kernels
  virtual void Update_(State& S) override;

 protected:
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override
  {
    AMANZI_ASSERT(false);
  }

 protected:
  Key elev_key_;
  Key pd_key_;
  Teuchos::RCP<Functions::MeshFunction> func_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorBCSeepageFaceHead> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
