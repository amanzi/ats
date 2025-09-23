/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! Macropore water retention model.
/*!

An exponential model basd on the equations of... (document me!)

.. math::
  ...


.. _WRM-macropore-spec
.. admonition:: WRM-macropore-spec

   * `"region`" ``[string]`` Region to which this applies
   * `"macropore exponent a [-]`" ``[double]`` The exponent alpha above.
   * `"residual saturation [-]`" ``[double]`` Limit at pc = inf
   * `"smoothing parameter [Pa]`" ``[double]`` **10** Exponential smoothing factor dP above.
   * `"water entry pressure [Pa]`" ``[double]`` Entry pressure Pwe above.

 */


#ifndef ATS_FLOWRELATIONS_WRM_MACROPORE_
#define ATS_FLOWRELATIONS_WRM_MACROPORE_

#include "Teuchos_ParameterList.hpp"
#include "Spline.hh"

#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class WRMMacropore : public WRM {
 public:
  explicit WRMMacropore(Teuchos::ParameterList& plist);

  // required methods from the base class
  double k_relative(double saturation);
  double d_k_relative(double saturation);
  double saturation(double pc);
  double d_saturation(double pc);
  double capillaryPressure(double saturation);
  double d_capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }
  //double suction_head(double saturation);
  //double d_suction_head(double saturation);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList& plist_;

  double a_; // Macropore parameters: a, dP, Pwe
  double dP_;
  double Pwe_;
  double sr_; // Macropore residual saturation

  // int function_;
  // double s0_;  // regularization threshold in saturation
  // Amanzi::Utils::Spline fit_kr_;

  // double pc0_;
  // Amanzi::Utils::Spline fit_s_;


  static Utils::RegisteredFactory<WRM, WRMMacropore> factory_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
