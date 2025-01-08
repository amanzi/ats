/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  ATS

  Saturated Vapor Pressure for vapor over water or ice, Sonntag (1990)

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
