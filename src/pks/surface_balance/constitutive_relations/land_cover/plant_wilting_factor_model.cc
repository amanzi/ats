/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Plant wilting factor provides a moisture availability-based limiter on transpiration.
#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "plant_wilting_factor_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
PlantWiltingFactorModel::PlantWiltingFactorModel(const LandCover& lc) : lc_(lc) {}

// main method
double
PlantWiltingFactorModel::PlantWiltingFactor(double pc) const
{
  return lc_.stomata_closed_water_potential < pc ?
           0. :
           (pc < lc_.stomata_open_water_potential ?
              1. :
              ((-pc + lc_.stomata_closed_water_potential) /
               (lc_.stomata_closed_water_potential - lc_.stomata_open_water_potential)));
}

double
PlantWiltingFactorModel::DPlantWiltingFactorDCapillaryPressureGasLiq(double pc) const
{
  return lc_.stomata_closed_water_potential < pc ?
           0. :
           (pc < lc_.stomata_open_water_potential ?
              0. :
              (-1 / (lc_.stomata_closed_water_potential - lc_.stomata_open_water_potential)));
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
