/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! WRMClappHornberger : water retention model using Clapp-Hornberger's parameterization

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov), F.-M. Yuan (yuanf@ornl.gov)
*/

/*!

Clapp-Hornberger's water retention curve.
NOTES:
(1) currently not yet include the parabolic section (inflection point - saturated)
    in order to be consistent with ELM's embedded Clapp-Hornberger algorithm, by default.
    IF 'smoothing' option on, it will be available for both near-saturation and extreme-dry saturation.
(2) Clapp-Hornberger WRM are very similar as Brook-Corey model.
    Clapp-Hornberger:    pc   = smpsat * S**(-bsw)
                         Krel = S**(3+2*bsw)

    Brook-Corey:         pc   =  (Se**(-1.0/lambda))/alpha, with Se=(S-Sr)/(1-Sr)
                         Krel = Se**(3+2/lamda), when soil wetting

    where, pc is capillary potential (Pa), Krel is relative permessivity, S is liq. water saturation, Sr is residual S, and Se effective S.
           AND,
           smpsat is soil matrical potential at saturated (positive, Clapp-Hornberger 'PSIs'), bsw is Clapp-Hornberger 'b';
           alpha and lamda are two parameters for Brook-Corey model.

    Assuming: Sr = 0.0, then: alpha = 1.0/smpsat, lamda=1.0/bsw

 (3) variable numerical ranges in this wrm. Dry-end: s<sx_; Near-saturated: s>s0_
   // s  ranges: -- 0 ---- sr ---- sx_ ----- s0_ -- (1.0-sr_) -- 1.0 --
   // se       : --------- 0. ----------------------- 1.0     ---------
  //  pc       : -------- Inf ---- pcx_ ---- pc0_ --- 0.0     ---------
  //  krel     : ----------------- 0.0  ------------- 1.0     ---------


.. _WRM-Clapp-Hornberger-spec
.. admonition:: WRM-Clapp-Hornberger-spec

    * `"region`" ``[string]`` Region to which this applies
    * `"Clapp Hornberger bsw [-]`" ``[double]`` Clapp-Hornberger's b parameter
    * `"Clapp Hornberger smpsat [Pa]`" ``[double]`` Clapp-Hornberger's minimum soil suction at saturated (air entry suction)
    * `"residual saturation [-]`" ``[double]`` **0.0** (if 0.0, it's exactly in Clapp-Hornberger Equations)
    * `"near-saturation inflection point interval [saturation]`" ``[double]`` **0.08** (inflection point, Si, in C-H 1978 paper. if > 0.0, a hyperbola eq. will be applied on near-saturation end. By default 0.08 and can be off by 0)
    * `"near-saturation smoothing starting point [Pa]`" ``[double]`` **0.0** (inflection point in capillary pressure and applying derivative-smoothing approach for until pc=0. So only if Si above not set and should be greater than 'smpsat'. If neither on, pure CH wrm used)
    * `"dry-end smoothing starting point [Pa]`" ``[double]`` **1.0e10** (if PSi over this point, PCx, either truncation or derivative-smoothing will be applied on WR curve dry-end. By default 1.e10 Pa and can be off by 0)
    * `"Krel function name`" ``[string]`` `"`" (hard-weired, i.e. user cannot modify it).

Example:

.. code-block:: xml

    <ParameterList name="soilA" type="ParameterList">
      <Parameter name="region" type="string" value="horizonA" />
      <Parameter name="WRM type" type="string" value="Clapp Hornberger" />
      <Parameter name="Clapp Hornberger smpsat [Pa]" type="double" value="4689.0" />
      <Parameter name="Clapp Hornberger bsw [-]" type="double" value="5.39" />
      <Parameter name="residual saturation [-]" type="double" value="0.0"/>
      <Parameter name="near-saturation injection point interval [saturation]" type="double" value="0.08" />
      <Parameter name="near-saturation smoothing starting point [Pa]" type="double" value="0.0" />
      <Parameter name="dry-end smoothing starting point [Pa]" type="double" value="1.0e10" />
    </ParameterList>

*/

#ifndef ATS_FLOWRELATIONS_WRM_CLAPP_HORNBERGER_
#define ATS_FLOWRELATIONS_WRM_CLAPP_HORNBERGER_

#include "Teuchos_ParameterList.hpp"
#include "Spline.hh"

#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMClappHornberger : public WRM {

public:
  explicit WRMClappHornberger(Teuchos::ParameterList& plist);

  // required methods from the base class
  double k_relative(double saturation);
  double d_k_relative(double saturation);
  double saturation(double pc);
  double d_saturation(double pc);
  double capillaryPressure(double saturation);
  double d_capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }
  double suction_head(double saturation);
  double d_suction_head(double saturation);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList& plist_;

  // Clapp-Hornberger parameters: smpsat, bsw
  double smpsat_;
  double bsw_;
  double sr_ = 0.0;  // residual saturation, by default is 0.0
  double m_;
  double n_;

  double s0_;       // near-saturation inflection point
  double pc0_;
  bool wet_step_smoothing_;   // sigmoid-smoothed step (Heaveside) function: (H1-H0)/[1+exp(-k*(x-(x1+x0)/2))]
  double sigmoid_one_eps_;     // allowable episilon for H=1 (i.e. gap from H=1)
  double sigmoid_k0_;         // exponential coefficient K
  double sigmoid_half0_;      // mid-point of x1~x0, at which H=1 and H=0 respectively
  Amanzi::Utils::Spline fit_kr0_;

  double pcx_;      // dry-end smoothing starting point
  double sx_;
  bool dry_step_smoothing_;
  double sigmoid_k1_;
  double sigmoid_half1_;
  //Amanzi::Utils::Spline fit_sx_;



  static Utils::RegisteredFactory<WRM,WRMClappHornberger> factory_;
};

} //namespace
} //namespace

#endif
