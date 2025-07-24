/*
  Copyright 2010–202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Saubhagya Rathore (rathoress@ornl.gov)

*/

//! Simulates water movement through culverts by instant transfer of water from inlet to outlet region.
/*!

Flow is calculated using standard culvert hydraulics, considering both inlet-controlled and outlet-controlled regimes.

Implements the following culvert hydraulics equations:
- **Inlet control**:

  .. math::
     Q_{\text{inlet}} = N_b \, C \, A \, \sqrt{2 g h_i}

  where:

  - :math:`N_b` = number of barrels *(unitless)*
  - :math:`C` = discharge coefficient *(unitless)*
  - :math:`A` = culvert cross-sectional area *(m^2)*
  - :math:`h_i` = head at culvert inlet *(m)*
  - :math:`g` = acceleration due to gravity *(9.81 m/s^2)*
  - :math:`Q_{\text{inlet}}` = inlet discharge *(m^3/s)*

- **Outlet control**:

  .. math::
     Q_{\text{outlet}} = N_b \, C \, A \, \sqrt{ \frac{2 g h_o}{k} }

  where:

  - :math:`h_o` = head difference between inlet and outlet *(m)*
  - :math:`k = 1.5 + \frac{29 n^2 L}{R^{4/3}}` *(unitless resistance term)*
  - :math:`n` = Manning's roughness coefficient *(s/m^{1/3})*
  - :math:`L` = culvert length *(m)*
  - :math:`R` = hydraulic radius *(m)*
  - :math:`Q_{\text{outlet}}` = outlet discharge *(m^3/s)*

- **Blended total discharge**:

  .. math::
     Q = \frac{Q_{\text{inlet}} \, Q_{\text{outlet}}}{\sqrt{Q_{\text{inlet}}^2 + Q_{\text{outlet}}^2 + \epsilon}}

  where:

  - :math:`\epsilon` = small positive constant to avoid divide-by-zero
  - :math:`Q` = total culvert discharge *(m^3/s)*

The resulting :math:`Q` is used to compute area-weighted water removal at the inlet and volume-weighted water addition at the outlet.


`"evaluator type`" = `"culvert`"

.. _evaluator-culvert-spec:
.. admonition:: evaluator-culvert-spec

   * `"culvert inlet region"`" ``[str]`` Region of cells where culvert flow is taken out.
   * `"culvert outlet region"`" ``[str]`` Region of cells where culvert flow is introduced.
   * `"number of barrels"`" ``[int]`` Number of culvert barrels, default is 1.
   * `"culvert length"`" ``[double]`` Length of the culvert in meters, default is 10.
   * `"culvert diameter"`" ``[double]`` Diameter of the culvert in meters, default is 1.
   * `"culvert roughness coefficient"`" ``[double]`` Manning's roughness coefficient for the culvert, default is 0.013.
   * `"culvert discharge coefficient"`" ``[double]`` Discharge coefficient for the culvert, default is 0.6.

   KEYS:
   - `"cell volume"`" **DOMAIN-cell_volume**
   - `"ponded depth"`" **DOMAIN-ponded_depth**
   - `"potential"`" **DOMAIN-pres_elev** (stage or water surface elevation)
   - `"elevation"`" **DOMAIN-elevation**
   - `"water content"`" **DOMAIN-water_content**
   - `"molar density liquid"`" **DOMAIN-molar_density_liquid**

Example:

.. code-block:: xml

      <ParameterList name="surface-culvert_flow" type="ParameterList">
        <Parameter name="evaluator type" type="string" value="culvert"/>
        <Parameter name="culvert inlet region" type="string" value="culvert inlet"/>
        <Parameter name="culvert outlet region" type="string" value="culvert outlet"/>
        <Parameter name="number of barrels" type="int" value="1"/>
        <Parameter name="culvert length" type="double" value="10.0"/>
        <Parameter name="culvert diameter" type="double" value="0.25"/>
        <Parameter name="culvert roughness coefficient" type="double" value="0.013"/>
        <Parameter name="culvert discharge coefficient" type="double" value="0.6"/>
      </ParameterList>

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "FunctionFactory.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class SurfCulvertEvaluator : public EvaluatorSecondaryMonotypeCV {
 public:
  explicit SurfCulvertEvaluator(Teuchos::ParameterList& plist);
  SurfCulvertEvaluator(const SurfCulvertEvaluator& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SurfCulvertEvaluator(*this));
  }

  // virtual void EnsureCompatibility(State& S) override;
  virtual bool IsDifferentiableWRT(const State& S,
                                   const Key& wrt_key,
                                   const Tag& wrt_tag) const override
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
                                          const std::vector<CompositeVector*>& result) override {};

 protected:
  Key domain_;
  Key cv_key_;
  Key pd_key_;
  Key pe_key_;
  Key elev_key_;
  Key liq_den_key_;
  Key wc_key_;

  std::string culvert_inlet_region_;
  std::string culvert_outlet_region_;
  int Nb_;
  double L_;
  double D_;
  double L_feet_;
  double D_feet_;
  double n_;
  double C_;
  bool allow_reverse_flow_;

 private:
  static Utils::RegisteredFactory<Evaluator, SurfCulvertEvaluator> reg_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
