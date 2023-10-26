/*
  Copyright 2010–202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Saubhagya Rathore (rathoress@ornl.gov)

*/

//! Simulates water movement through culverts by instant transfer of water from inlet to outlet region.
//! Flow is calculated using standard culvert hydraulics, considering both inlet-controlled and outlet-controlled regimes.
/*!
.. _evaluator-surface-culvert-spec:
.. admonition:: evaluator-surface-culvert-spec

   * `"culvert inlet region`" ``[str]`` Region of cells where culvert flow is taken out.
   * `"culvert outlet region`" ``[str]`` Region of cells where culvert flow is introduced.
   * `"number of barrels`" ``[int]`` Number of culvert barrels, default is 1.
   * `"culvert length`" ``[double]`` Length of the culvert in meters, default is 10.
   * `"culvert diameter`" ``[double]`` Diameter of the culvert in meters, default is 1.
   * `"culvert roughness coefficient`" ``[double]`` Manning's roughness coefficient for the culvert, default is 0.013.
   * `"culvert discharge coefficient`" ``[double]`` Discharge coefficient for the culvert, default is 0.6.

   KEYS:
   - `"cell volume`" **DOMAIN-cell_volume** 
   - `"ponded depth`" **DOMAIN-ponded_depth** 
   - `"potential`" **DOMAIN-pres_elev** (stage or water surface elevation)
   - `"elevation`" **DOMAIN-elevation** 
   - `"water content`" **DOMAIN-water_content** 
   - `"molar density liquid`" **DOMAIN-molar_density_liquid**

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

Implements the following culvert hydraulics equations:

   - Inlet control:
     \f[
     Q_{inlet} = N_b C A \sqrt{2g h_i}
     \f]
     where:
       - \( N_b \) = number of barrels  
       - \( C \) = discharge coefficient  
       - \( A \) = culvert cross-sectional area  
       - \( h_i \) = head at culvert inlet  
       - \( g \) = gravity

   - Outlet control:
     \f[
     Q_{outlet} = N_b C A \sqrt{ \frac{2g h_o}{k} }
     \f]
     where:
       - \( h_o \) = head difference between inlet and outlet  
       - \( k = 1.5 + \frac{29 n^2 L}{R^{4/3}} \) (Manning-based resistance term)  
       - \( n \) = Manning’s roughness  
       - \( L \) = culvert length  
       - \( R \) = hydraulic radius

   - Blended total discharge:
     \f[
     Q = \frac{Q_{inlet} Q_{outlet}}{\sqrt{Q_{inlet}^2 + Q_{outlet}^2 + \epsilon}}
     \f]
     where \( \epsilon \) is a small number to avoid divide-by-zero.

   The resulting \( Q \) is used to compute area-weighted water removal at the inlet and volume-weighted water addition at the outlet.
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
