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

#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class PlantWiltingFactorModel {
 public:
  explicit PlantWiltingFactorModel(const LandCover& lc);

  double PlantWiltingFactor(double pc) const;

  double DPlantWiltingFactorDCapillaryPressureGasLiq(double pc) const;

 protected:
  const LandCover& lc_;
};

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
