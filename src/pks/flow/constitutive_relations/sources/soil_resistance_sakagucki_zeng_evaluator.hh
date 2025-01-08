/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*!
Evaluates the soil resistance at top cells through the Sakagucki-Zeng model
referred to Sakagucki and Zeng (2009).

`"evaluator type`" = `"soil resistance, Sakagucki-Zeng`"

.. _soil_resistance_sakagucki_zeng_evaluator-spec
.. admonition:: soil_resistance_sakagucki_zeng_evaluator-spec

  * `"model parameters`" ``[string]`` **WRM parameters** ``[WRM-typedinline-spec-list]``
  Soil resistance based on Sakagucki-Zeng method uses soil properties defined in
  `"WRM parameters`" which is given through `"model parameters`" under state.

    - If `"van Genuchten`" is used for WRM, either `"van Genuchten n [-]`"
    or `"van Genuchten m [-]`" will be used to determine Clapp-Hornberger-b
    through method 2 in Ma et al. (1999). Originally this method is from
    Lenhard et al. (1989).

    - If `"Brooks-Corey`" is used for WRM, `"Brooks-Corey lambda [-]`" will
    be used to determine Clapp-Hornberger-b, which is the reciprocal of
    `"Brooks-Corey lambda [-]`".

    - `"residual saturation [-]`" ``[double]`` **0.0**

    - `"dessicated zone thickness [m]`" ``[double]`` **0.1**

  KEYS:

  - `"gas saturation`" of top cells
  - `"porosity`" of top cells

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
          <Parameter name="smoothing interval width [saturation]" type="double" value="0.05" />
          <Parameter name="dessicated zone thickness [m]" type="double" value="0.1" />
        </ParameterList>
      </ParameterList>
    </ParameterList>
    <ParameterList name="evaluators" type="ParameterList">
      <ParameterList name="surface-soil_resistance" type="ParameterList">
        <Parameter name="evaluator type" type="string" value="soil resistance, Sakagucki-Zeng" />
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
#include "soil_resistance_model_partition.hh"

namespace Amanzi {
namespace Flow {

class SoilResistanceSakaguckiZengEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SoilResistanceSakaguckiZengEvaluator(Teuchos::ParameterList& plist);
  SoilResistanceSakaguckiZengEvaluator(const SoilResistanceSakaguckiZengEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

  Teuchos::RCP<SoilResistanceModelPartition> get_Models() { return models_; }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  Key sat_gas_key_;
  Key poro_key_;

  Teuchos::RCP<SoilResistanceModelPartition> models_;

 private:
  static Utils::RegisteredFactory<Evaluator, SoilResistanceSakaguckiZengEvaluator> fac_;
};

} // namespace Flow
} // namespace Amanzi
