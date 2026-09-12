/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Plant wilting factor provides a moisture availability-based limiter on transpiration.
/*!

This implements the simple water-based limiter, given pressure in [Pa]:

.. math:
   \beta =  \frac{p_{closed} - p}{p_{closed} - p_{open}}

where p is the capillary pressure or water potential, and closed and open
indicate the values at which stomates are fully open or fully closed (the
wilting point).  These two parameters are provided by the LandCover object.

All parameters are in units of [Pa], and are positive (water potential).

*/

#pragma once

#include "Kokkos_Core.hpp"
#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class PlantWiltingFactorModel {
 public:
  KOKKOS_INLINE_FUNCTION
  explicit PlantWiltingFactorModel(const LandCover& lc)
    : pc_closed_(lc.stomata_closed_capillary_pressure),
      pc_open_(lc.stomata_open_capillary_pressure)
  {}

  KOKKOS_INLINE_FUNCTION
  double PlantWiltingFactor(double pc) const
  {
    return pc_closed_ < pc ? 0. :
                             (pc < pc_open_ ? 1. : ((-pc + pc_closed_) / (pc_closed_ - pc_open_)));
  }

  KOKKOS_INLINE_FUNCTION
  double DPlantWiltingFactorDCapillaryPressureGasLiq(double pc) const
  {
    return pc_closed_ < pc ? 0. : (pc < pc_open_ ? 0. : (-1 / (pc_closed_ - pc_open_)));
  }

 protected:
  double pc_closed_;
  double pc_open_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
