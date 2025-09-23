/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates Kr from surface into the subsurface

*/

#ifndef AMANZI_FLOWRELATIONS_SURFACE_KR_MODEL_
#define AMANZI_FLOWRELATIONS_SURFACE_KR_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class SurfaceRelPermModel {
 public:
  virtual ~SurfaceRelPermModel() {}
  virtual bool TemperatureDependent() = 0;
  virtual double SurfaceRelPerm(double uf, double h) = 0;
  virtual double DSurfaceRelPermDUnfrozenFraction(double uf, double h) = 0;
  virtual double DSurfaceRelPermDPondedDepth(double uf, double h) = 0;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
