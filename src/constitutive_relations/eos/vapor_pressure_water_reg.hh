/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  Saturated Vapor Pressure for vapor over water or ice, Sonntag (1990)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include <cmath>
#include "errors.hh"
#include "vapor_pressure_water.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<VaporPressureRelation, VaporPressureWater>
  VaporPressureWater::factory_("water vapor over water/ice");

} // namespace Relations
} // namespace Amanzi
