/*
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
PlantWiltingFactorModel::PlantWiltingFactorModel(const LandCover& lc)
  : lc_(lc) {}

// main method
double
PlantWiltingFactorModel::PlantWiltingFactor(double pc) const
{
  return lc_.stomata_closed_mafic_potential < pc ? 0. :
    (pc < lc_.stomata_open_mafic_potential ? 1. :
     ((-pc + lc_.stomata_closed_mafic_potential)/(lc_.stomata_closed_mafic_potential - lc_.stomata_open_mafic_potential)));
}

double
PlantWiltingFactorModel::DPlantWiltingFactorDCapillaryPressureGasLiq(double pc) const
{
  return lc_.stomata_closed_mafic_potential < pc ? 0. :
    (pc < lc_.stomata_open_mafic_potential ? 0. :
     (-1/(lc_.stomata_closed_mafic_potential - lc_.stomata_open_mafic_potential)));
}

} //namespace
} //namespace
} //namespace
