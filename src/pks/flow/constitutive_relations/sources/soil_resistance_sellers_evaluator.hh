/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*!
Evaluates the soil resistance at top cells through the Sellers model
referred to Sellers et al. (1992).

.. math::

      R_{soil} = \mathrm{exp}(8.206 - 4.255 s_l)


`"evaluator type`" = `"soil resistance, Sellers`"

.. _evaluator-soil-resistance-sellers-spec
.. admonition:: evaluator-soil-resistance-sellers-spec

   KEYS:

   - `"liquid saturation`" of top cells

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
        <Parameter name="evaluator type" type="string" value="soil resistance, Sellers" />
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

namespace Amanzi {
namespace Flow {

class SoilResistanceSellersModel;

class SoilResistanceSellersEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SoilResistanceSellersEvaluator(Teuchos::ParameterList& plist);
  SoilResistanceSellersEvaluator(const SoilResistanceSellersEvaluator& other) = default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;
  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& result) override;

  virtual void EnsureCompatibility_ToDeps_(State& S) override;

 protected:
  Key sat_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, SoilResistanceSellersEvaluator> fac_;
};

} // namespace Flow
} // namespace Amanzi
