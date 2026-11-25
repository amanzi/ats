/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
Evaluates actual bare-soil evaporation from potential evaporation using the
Lee-Pielke downregulation, a smooth cosine ramp of the top grid cell's
liquid saturation between full evaporation and a water-retention-model-based
field capacity.

.. math::

   E = E_{potential} \, \beta \quad \text{if } E_{potential} > 0, \text{ else } E_{potential}

   \beta = \begin{cases}
     1 & (1 - s_g) \geq s_{fc} \\
     0.25 \left[ 1 - \cos\left( \pi \frac{1 - s_g}{s_{fc}} \right) \right]^2 & \text{otherwise}
   \end{cases}

where :math:`s_g` is the top subsurface cell's gas saturation and
:math:`s_{fc}` is the liquid saturation at field capacity for that column's
water retention model (see ``WRM::field_capacity()``).

`"evaluator type`" = `"evaporation downregulation, Lee-Pielke`"

.. _evaluator-evaporation-downregulation-lee-pielke-spec:
.. admonition:: evaluator-evaporation-downregulation-lee-pielke-spec

   * `"model parameters`" ``[string]`` **WRM parameters** Field capacity is
     computed from the water retention models defined in `"WRM parameters`"
     under state (see `"field capacity relative permeability [-]`" on each
     WRM).

   KEYS:

   - `"potential evaporation`"
   - `"gas saturation`" Note only the column's top grid cell is used

Example:

.. code-block:: xml

  <ParameterList name="state" type="ParameterList">
    <ParameterList name="model parameters" type="ParameterList">
      <ParameterList name="WRM parameters" type="ParameterList">
        <ParameterList name="domain" type="ParameterList">
          <Parameter name="region" type="string" value="domain" />
          <Parameter name="wrm type" type="string" value="van Genuchten" />
          <Parameter name="van Genuchten alpha [Pa^-1]" type="double" value="2e-05" />
          <Parameter name="van Genuchten n [-]" type="double" value="1.58" />
          <Parameter name="residual saturation [-]" type="double" value="0.2" />
          <Parameter name="field capacity relative permeability [-]" type="double" value="1.e-3" />
        </ParameterList>
      </ParameterList>
    </ParameterList>
    <ParameterList name="evaluators" type="ParameterList">
      <ParameterList name="surface-evaporation" type="ParameterList">
        <Parameter name="evaluator type" type="string" value="evaporation downregulation, Lee-Pielke" />
        <Parameter name="model parameters" type="string" value="WRM parameters" />
      </ParameterList>
      ...
    </ParameterList>
    ...
  </ParameterList>

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "wrm_partition.hh"

namespace Amanzi {
namespace Flow {

class EvaporationDownregulationLeePielkeEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit EvaporationDownregulationLeePielkeEvaluator(Teuchos::ParameterList& plist);
  EvaporationDownregulationLeePielkeEvaluator(
    const EvaporationDownregulationLeePielkeEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<WRMPartition> get_WRMs() { return wrms_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  Key pot_evap_key_;
  Key sat_gas_key_;

  Teuchos::RCP<WRMPartition> wrms_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaporationDownregulationLeePielkeEvaluator> fac_;
};

} // namespace Flow
} // namespace Amanzi
