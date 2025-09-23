/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Linear interpolant of thermal conductivity.
------------------------------------------------------------------------- */

#include <cmath>
#include "dbc.hh"
#include "thermal_conductivity_threephase_volume_averaged.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

ThermalConductivityThreePhaseVolumeAveraged::ThermalConductivityThreePhaseVolumeAveraged(
  Teuchos::ParameterList& plist)
  : plist_(plist)
{
  InitializeFromPlist_();
};

double
ThermalConductivityThreePhaseVolumeAveraged::ThermalConductivity(double poro,
                                                                 double sat_liq,
                                                                 double sat_ice,
                                                                 double temp)
{
  AMANZI_ASSERT(std::abs(1 - sat_liq - sat_ice) < 1.e-10);
  return (1 - poro) * k_soil_ + poro * sat_liq * k_liquid_ + poro * sat_ice * k_ice_ +
         poro * (1 - sat_liq - sat_ice) * k_gas_;
};

void
ThermalConductivityThreePhaseVolumeAveraged::InitializeFromPlist_()
{
  k_soil_ = plist_.get<double>("thermal conductivity of soil [W m^-1 K^-1]");
  k_ice_ = plist_.get<double>("thermal conductivity of ice [W m^-1 K^-1]");
  k_liquid_ = plist_.get<double>("thermal conductivity of liquid [W m^-1 K^-1]");
  k_gas_ = plist_.get<double>("thermal conductivity of gas [W m^-1 K^-1]");
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
