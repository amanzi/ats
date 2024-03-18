/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

//! WRMBrooksCorey : water retention model using Brooks-Corey's parameterization
/*!

Brooks-Corey's water retention curve, typically used to determine relative
permeability under freezing conditions by converting van Genuchten parameters
to Brooks-Corey parameters.

.. _WRM-Brooks-Corey-spec
.. admonition:: WRM-Brooks-Corey-spec

   * `"region`" ``[string]`` Region to which this applies
   * `"Brooks-Corey lambda [-]`" ``[double]`` Brooks-Corey's lambda
   * `"Brooks-Corey saturated matric suction [Pa]`" ``[double]`` Brooks-Corey
     saturated matric suction in Pa
   * `"residual saturation [-]`" ``[double]`` **0.0**
   * `"smoothing interval width [saturation]`" ``[double]`` **0.0**

Example:

.. code-block:: xml

    <ParameterList name="moss" type="ParameterList">
      <Parameter name="region" type="string" value="moss" />
      <Parameter name="WRM type" type="string" value="Brooks-Corey" />
      <Parameter name="Brooks-Corey lambda [-]" type="double" value="0.5" />
      <Parameter name="Brooks-Corey saturated matric suction [Pa]" type="double" value="1.e3" />
      <Parameter name="residual saturation [-]" type="double" value="0.0" />
      <Parameter name="smoothing interval width [saturation]" type="double" value=".05" />
    </ParameterList>

*/

#ifndef ATS_FLOWRELATIONS_WRM_BROOKS_COREY_
#define ATS_FLOWRELATIONS_WRM_BROOKS_COREY_

#include "Teuchos_ParameterList.hpp"
#include "Spline.hh"

#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMBrooksCorey : public WRM {
 public:
  explicit WRMBrooksCorey(Teuchos::ParameterList& plist);

  // required methods from the base class
  double k_relative(double saturation);
  double d_k_relative(double saturation);
  double saturation(double pc);
  double d_saturation(double pc);
  double capillaryPressure(double saturation);
  double d_capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList& plist_;

  double lambda_, p_sat_; // Brooks and Corey parameters: lambda, p_sat
  double sr_;             // residual saturation
  double b_;              // frequently used constant, clapp-hornberger b = 1/lambda

  double s0_; // regularization threshold in saturation
  Amanzi::Utils::Spline fit_kr_;

  static Utils::RegisteredFactory<WRM, WRMBrooksCorey> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
