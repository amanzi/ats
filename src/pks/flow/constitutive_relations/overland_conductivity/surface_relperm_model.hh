
/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates Kr from surface into the subsurface

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_SURFACE_KR_MODEL_
#define AMANZI_FLOWRELATIONS_SURFACE_KR_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
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
} // namespace Amanzi

#endif
