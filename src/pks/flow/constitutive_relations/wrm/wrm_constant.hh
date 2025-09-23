/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! Constant water retention model.
/*!

A constant-valued WRM.

.. _WRM-constant-spec
.. admonition:: WRM-constant-spec

    * `"region`" ``[string]`` Region to which this applies
    * `"constant value [-]`" ``[double]`` **1.0** Constant relative permeability value.
    * `"residual saturation [-]`" ``[double]`` **1.0**  Constant saturation value.

Example:

.. code-block:: xml

    <ParameterList name="moss" type="ParameterList">
      <Parameter name="region" type="string" value="moss" />
      <Parameter name="WRM type" type="string" value="constant" />
    </ParameterList>


 */

#ifndef ATS_FLOWRELATIONS_WRM_CONSTANT_
#define ATS_FLOWRELATIONS_WRM_CONSTANT_

#include "Teuchos_ParameterList.hpp"
#include "Spline.hh"

#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class WRMConstant : public WRM {
 public:
  explicit WRMConstant(Teuchos::ParameterList& plist);

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

  double a_;
  double sr_;

  static Utils::RegisteredFactory<WRM, WRMConstant> factory_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
