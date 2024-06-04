/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
A boundary condition evaluator for seepage faces.
*/

/*!

`"evaluator type`" == `"flow bcs seepage`"

.. _flow-bcs-seepage-spec:
.. admonition:: flow-bcs-seepage-spec

   * `"constant in time`" ``[bool]`` **false** If true, only evaluate the
     functions once as they are time-independent.
   * `"pressure`" ``[mesh-function-spec-list]`` The pressure above which water can seep.
   * `"water flux`" ``[mesh-function-spec-list]`` The (outward) flow rate when
     pressure is below the seepage pressure.  Typically this is 0, but is
     sometimes negative for infiltration when not seeping.

*/

#pragma once

#include "MeshFunction.hh"
#include "EvaluatorSecondary.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class EvaluatorFlowBCsSeepage : public EvaluatorSecondary {
 public:
  explicit EvaluatorFlowBCsSeepage(const Teuchos::RCP<Teuchos::ParameterList>& plist);
  EvaluatorFlowBCsSeepage(const EvaluatorFlowBCsSeepage& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual Evaluator& operator=(const Evaluator& other) override;

  EvaluatorFlowBCsSeepage& operator=(const EvaluatorFlowBCsSeepage& other);

  virtual std::string getType() const override { return "independent variable patch"; }

  virtual void EnsureCompatibility(State& S) override;

 protected:
  virtual void Update_(State& S) override;

 protected:
  Teuchos::RCP<Functions::MeshFunction> func_;
  std::string function_outer_name_, function_inner_name_;
  AmanziMesh::Entity_kind entity_kind_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorFlowBCsSeepage> fac_;
};

} // namespace Amanzi
