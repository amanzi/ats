/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Bo Gao (gaob@ornl.gov)
*/

/*!
Evaluates the slope stability as a factor of safety for subsurface cells using
TRIGRS, referred to Baum et al. (2008).

.. math::

   FS = \frac{\tan \varphi'}{\tan \sigma}
      + \frac{c' - \chi (p - p_{atm}) \tan \varphi'}{\rho_{bulk} g d \sin \sigma \cos \sigma}

   \chi = \frac{s_l - s_r}{1 - s_r}
   \rho_{bulk} = \rho_s (1 - \phi) + \rho_{liq} \phi * s_l

- :math:`\varphi'` effective internal friction angle ``[rad]``
- :math:`\sigma` slope angle ``[rad]``
- :math:`c'` effective soil cohesion ``[Pa]``
- :math:`p` pressure ``[Pa]``, :math:`p_{atm}` atmospheric pressure ``[Pa]``
- :math:`s_l` liquid saturation ``[-]``, :math:`s_r` residual saturation ``[-]``
- :math:`\rho_s` soil grain density ``[kg m^-3]``
- :math:`\phi` porosity ``[-]``
- :math:`\rho_{liq}` liquid density ``[kg m^-3]``
- :math:`d` depth of cells ``[m]``


`"evaluator type`" = `"slope stability, TRIGRS`"

.. _evaluator-slope-stability-trigrs-spec:
.. admonition:: evaluator-slope-stability-trigrs-spec

   * `"model parameters`" ``[string]`` **WRM parameters** Slope stability based
     on TRIGRS method uses soil properties defined in `"WRM parameters`"
     which is given through `"model parameters`" under state.

     - `"residual saturation [-]`" ``[double]`` **0.0**
     - `"minimum slope [rad]`" ``[double]`` **1e-6** Cells with a slope angle
       smaller than this are treated as flat.
     - `"maximum slope stability [-]`" ``[double]`` **1e10** Factor of safety
       value assigned to flat cells.

   KEYS:

   - `"effective internal friction`" **DOMAIN-effective_internal_friction** :math:`\phi'`
   - `"effective soil cohesion`" **DOMAIN-effective_soil_cohesion** :math:`c'`
   - `"pressure`"
   - `"saturation`"
   - `"depth`"
   - `"density rock`"
   - `"mass density liquid`"
   - `"porosity`"
   - `"slope`" **SURFACE_DOMAIN-slope_magnitude**

Example:

.. code-block:: xml

  <ParameterList name="state" type="ParameterList">
    <ParameterList name="model parameters" type="ParameterList">
      <ParameterList name="WRM parameters" type="ParameterList">
        <ParameterList name="domain" type="ParameterList">
          <Parameter name="region" type="string" value="domain" />
          ... ...
          <Parameter name="residual saturation [-]" type="double" value="0.2" />
          <Parameter name="minimum slope [rad]" type="double" value="1.e-6" />
          <Parameter name="maximum slope stability [-]" type="double" value="1.e10" />
        </ParameterList>
      </ParameterList>
    </ParameterList>
    <ParameterList name="evaluators" type="ParameterList">
      <ParameterList name="slope_stability" type="ParameterList">
        <Parameter name="evaluator type" type="string" value="slope stability, TRIGRS" />
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
#include "region_model_partition.hh"

namespace Amanzi {
namespace Flow {

class SlopeStabilityTrigrsModel;
using SlopeStabilityTrigrsModelPartition = RegionModelPartition<SlopeStabilityTrigrsModel>;

class SlopeStabilityTrigrsEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SlopeStabilityTrigrsEvaluator(Teuchos::ParameterList& plist);
  SlopeStabilityTrigrsEvaluator(const SlopeStabilityTrigrsEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  void InitializeFromPlist_();
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  Key domain_sub_, domain_surf_;
  Key internal_friction_key_;
  Key soil_cohesion_key_;
  Key pressure_key_;
  Key sat_key_;
  Key depth_key_;
  Key rho_r_key_;
  Key poro_key_;
  Key rhow_key_;
  Key slope_key_;

  Teuchos::RCP<SlopeStabilityTrigrsModelPartition> models_;

 private:
  static Utils::RegisteredFactory<Evaluator, SlopeStabilityTrigrsEvaluator> fac_;
};

} // namespace Flow
} // namespace Amanzi
