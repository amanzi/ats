/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Helper functions for the subgrid topography model. Note these evaluate equation 7 of Jan et al WRR 2018, and its derivative.
#pragma once

#include <cmath>

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Microtopography {

inline double
volumetricDepth(double depth, double del_max, double del_ex)
{
  if (depth == 0.) {
    return 0.;
  } else if (depth > del_max) {
    return depth - del_ex;
  } else {
    return std::pow(depth / del_max, 2) * (2 * del_max - 3 * del_ex) +
           std::pow(depth / del_max, 3) * (2 * del_ex - del_max);
  }
}

inline double
dVolumetricDepth_dDepth(double depth, double del_max, double del_ex)
{
  if (depth == 0.) {
    return 0.;
  } else if (depth > del_max) {
    return 1.;
  } else {
    return 2 * depth / del_max * (2 * del_max - 3 * del_ex) / del_max +
           3 * std::pow(depth / del_max, 2) * (2 * del_ex - del_max) / del_max;
  }
}

inline bool
validParameters(double del_max, double del_ex)
{
  return (2 * del_max >= 3 * del_ex) && (2 * del_ex >= del_max);
}


} // namespace Microtopography
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
