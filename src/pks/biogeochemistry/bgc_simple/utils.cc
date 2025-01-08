/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Chonggang Xu (cxu@lanl.gov)
*/

/*

Functor for evaluating QSat

*/

#include <cmath>
#include "utils.hh"

namespace Amanzi {
namespace BGC {

// returns the depth of the face whose cell above is thawed but cell
// below is frozen
double
PermafrostDepth(const Epetra_SerialDenseVector& SoilTArr,
                const Epetra_SerialDenseVector& SoilThicknessArr,
                double freeze_temp)
{
  int i = 0;
  double depth = 0.;
  int nSoilLayers = SoilTArr.Length();
  while (i < nSoilLayers && SoilTArr[i] > freeze_temp) {
    depth += SoilThicknessArr[i];
    i++;
  }
  return depth;
}

// returns the index of the top-most frozen cell, or ncells if all are
// frozen
int
PermafrostDepthIndex(const Epetra_SerialDenseVector& SoilTArr, double freeze_temp)
{
  int i = 0;
  int nSoilLayers = SoilTArr.Length();
  while (i < nSoilLayers && SoilTArr[i] > freeze_temp) i++;
  return i;
}

// This function calculate the effect of temperature on biological process.
double
TEffectsQ10(double Q10, double T, double refT)
{
  return std::pow(Q10, 0.1 * (T - refT));
}

} // namespace BGC
} // namespace Amanzi
