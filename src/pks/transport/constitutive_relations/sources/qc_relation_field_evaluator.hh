/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Phong V.V. Le (lepv@ornl.gov)
*/

//! Evaluates transport source (mass) from a field discharge
/*!

Mass sources into stream/river from a field discharge

.. _qc_relation_field_evaluator-spec:
.. admonition:: qc_relation_field_evaluator-spec

   `"function`" ``[function-spec]`` Function describing the relationship between field discharge (e.g. tile, groundwater) and solute mass running into river/stream

   KEYS:
   - `"cell volume`" **DOMAIN-cell_volume**
   - `"molar density liquid`" **DOMAIN-molar_density_liquid**
   - `"field source`" **DOMAIN-field_source** source
   - `"extensive`" ``[bool]`` checks if source is extensive. Default value is *false*.

Example

.. code-block:: xml

  <ParameterList name="RIVER_DOMAIN-field_sources">
    <Parameter name="evaluator type" type="string" value="q-c field" />
    <Parameter name="field source key" type="string" value="RIVER_DOMAIN-water_source_field" />
    <Parameter name="extensive" type="bool" value="false" />
    <ParameterList name="function" type="ParameterList">
      <ParameterList name="function-tabular" type="ParameterList">
        <Parameter name="x values" type="Array(double)" value="{0, 0.1, 1, 10}" />
        <Parameter name="y values" type="Array(double)" value="{0, 0.08, 0.1, 0.2}" />
        <Parameter name="forms" type="Array(string)" value="{linear, linear, linear}" />
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

class QCRelationFieldEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit QCRelationFieldEvaluator(Teuchos::ParameterList& plist);
  QCRelationFieldEvaluator(const QCRelationFieldEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new QCRelationFieldEvaluator(*this));
  }

  // virtual void EnsureCompatibility(State& S) override;
  virtual bool IsDifferentiableWRT(const State& S,
                                   const Key& wrt_key,
                                   const Tag& wrt_tag) const override
  {
    return false;
  }

 protected:
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override {};

 protected:
  Key domain_;
  Key cv_key_;
  Key molar_density_key_;
  Key tcc_key_;
  Key field_src_key_;
  bool extensive_;
  Teuchos::RCP<Function> QC_curve_;

 private:
  static Utils::RegisteredFactory<Evaluator, QCRelationFieldEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
