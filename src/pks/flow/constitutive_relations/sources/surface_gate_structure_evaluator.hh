/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Saubhagya Rathore (rathoress@ornl.gov)
           Ethan Coon (coonet@ornl.gov)
*/

//! Simulates gravity-driven water movement/diversions through gate structure.
/*!

This evaluator models gate structure inside 2D flow area. The gate curve is kept general and can be provided by the user.
Gate structure can be used to move water two canals, two storage areas or canal to storage area.

`"evaluator type`" = `"gate structure`"

.. _evaluator-surface-gate-structure-spec:
.. admonition:: evaluator-surface-gate-structure-spec

   * `"gate intake region`" ``[str]`` Region of cells where gate flow is taken out.
   * `"storage area region`" ``[str]`` Region of cells where gate flow is introduced.
   * `"gate stage close`" ``[double]`` The water surface elevation in the storage area that should trigger the gate close.
   * `"is ponded depth function`" ``[bool]`` If true, the gate efficiency curve is a function of ponded depth in the intake region rather than stage.
   * `"function" ``[function-tabular]`` This is a function/table of head and flow. Here, head is the stage or ponded depth in the intake region (upstream) of the gate structure.

   KEYS:
   - `"cell volume`" **DOMAIN-cell_volume** 
   - `"ponded depth`" **DOMAIN-ponded_depth** 
   - `"potential`" **DOMAIN-pres_elev** stage or water surface elevation
   - `"elevation`" **DOMAIN-elevation** 
   - `"water content`" **DOMAIN-water_content** 
   - `"molar liquid density`" **DOMAIN-molar_density_liquid** 

Example:

.. code-block:: xml

      <ParameterList name="surface-gate_flow" type="ParameterList">
        <Parameter name="evaluator type" type="string" value="gate structure"/>
        <Parameter name="gate intake region" type="string" value="gate intake area"/>
        <Parameter name="storage area region" type="string" value="detention pond"/>
        <Parameter name="gate close stage" type="double" value="7.50"/>
        <Parameter name="is ponded depth function" type="bool" value="true"/>
        <ParameterList name="function" type="ParameterList">
          <ParameterList name="function-tabular" type="ParameterList">
            <Parameter name="x values" type="Array(double)" value="{0.0, 0.18, 0.1840954, 0.20202941, 0.21996341, 0.23789742, 0.25583143}"/>
            <Parameter name="y values" type="Array(double)" value="{0.0, 0.0, 0.01542766, 0.07531978, 0.1352119, 0.19510401, 0.25499613}"/>
            <Parameter name="forms" type="Array(string)" value="{linear, linear, linear, linear, linear, linear}"/>
          </ParameterList>
        </ParameterList>
      </ParameterList>
*/


#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class SurfGateEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SurfGateEvaluator(Teuchos::ParameterList& plist);
  SurfGateEvaluator(const SurfGateEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SurfGateEvaluator(*this));
  }

  // virtual void EnsureCompatibility(State& S) override;
  virtual bool
  IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    // calculate of derivatives of this is a tricky thing to do, with
    // non-cell-local terms due to rescaling.  Just turn off derivatives
    // instead.
    return false;
  }
  
 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                        const Key& wrt_key,
                                        const Tag& wrt_tag,
                                        const std::vector<CompositeVector*>& result) override{};

 protected:
  Key domain_;
  Key cv_key_;
  Key pd_key_;
  Key pe_key_;
  Key elev_key_;
  Key liq_den_key_;
  Key wc_key_;

  std::string storage_region_;
  std::string gate_intake_region_;
  double stage_close_;
  bool is_ponded_depth_function_; 

  Teuchos::RCP<Function> Q_gate_;

 private:
  static Utils::RegisteredFactory<Evaluator, SurfGateEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

