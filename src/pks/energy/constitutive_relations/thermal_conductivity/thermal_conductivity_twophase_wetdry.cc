/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS


Simple model of two-phase thermal conductivity, based upon:

- Interpolation between saturated and dry conductivities via a Kersten number.
- Power-law Kersten number.

See ATS process model documentation's permafrost model for details.
------------------------------------------------------------------------- */

#include <cmath>
#include "thermal_conductivity_twophase_wetdry.hh"

namespace Amanzi {
namespace Energy {

// constructor
ThermalConductivityTwoPhaseWetDry::ThermalConductivityTwoPhaseWetDry(Teuchos::ParameterList& plist)
  : plist_(plist)
{
  InitializeFromPlist_();
};

// do the physics
double
ThermalConductivityTwoPhaseWetDry::ThermalConductivity(double poro, double sat_liq)
{
  double kersten = pow(sat_liq + eps_, alpha_);
  return k_dry_ + (k_wet_ - k_dry_) * kersten;
};

// initialization
void
ThermalConductivityTwoPhaseWetDry::InitializeFromPlist_()
{
  eps_ = plist_.get<double>("epsilon [-]", 1.e-10);
  alpha_ = plist_.get<double>("unsaturated alpha [-]", 1.0);
  k_dry_ = plist_.get<double>("thermal conductivity, dry [W m^-1 K^-1]");
  k_wet_ = plist_.get<double>("thermal conductivity, wet [W m^-1 K^-1]");
};

} // namespace Energy
} // namespace Amanzi
