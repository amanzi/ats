/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Data structures and parameters for calculating the snow / surface energy balance.
*/


#include <cmath>

#include "seb_physics_defs.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace SurfaceBalance {
namespace Relations {


Partition
Partitioner::CalcPartition(double ht_snow, double ht_pond, double unfrozen_fraction) const
{
  Partition part;

  if (ht_snow > snow_pen) {
    // incoming radiation does not penetrate the snow
    part.perSnow = 1.;
    part.perTundra = 0.;
    part.perWater = 0.;
    part.perIce = 0.;

  } else {
    part.perSnow = std::pow(ht_snow / snow_pen, 2);
    double remainder = 1. - part.perSnow;

    if (ht_pond > water_pen) {
      // incoming radiation does not penetrate snow and ponding
      part.perWater = unfrozen_fraction * remainder;
      part.perIce = remainder - part.perWater;
      part.perTundra = 0.;
    } else {
      double perPond = remainder * ht_pond / water_pen;
      part.perWater = unfrozen_fraction * perPond;
      part.perIce = (1. - unfrozen_fraction) * perPond;
      part.perTundra = remainder - perPond;
    }
  }
  return part;
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
