/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Chonggang Xu (cxu@lanl.gov)
*/

/*
Utility functions for Vegetation.

*/


#ifndef ATS_BGC_QSAT_HH_
#define ATS_BGC_QSAT_HH_

#include "Epetra_SerialDenseVector.h"

namespace Amanzi {
namespace ATS_Physics {
namespace BGC {


struct MetData {
  double qSWin;
  double tair;
  double windv;
  double wind_ref_ht;
  double vp_air;
  double CO2a;
  double lat;
};

double PermafrostDepth(const Epetra_SerialDenseVector& SoilTArr,
                       const Epetra_SerialDenseVector& SoilThicknessArr,
                       double freeze_temp);

int PermafrostDepthIndex(const Epetra_SerialDenseVector& SoilTArr, double freeze_temp);

// This function calculate the effect of temperature on biological process.
double TEffectsQ10(double Q10, double T, double refT);

} // namespace BGC
} // namespace ATS_Physics
} // namespace Amanzi


#endif
