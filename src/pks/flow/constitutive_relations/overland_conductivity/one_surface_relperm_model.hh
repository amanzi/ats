/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  No special limits as p_surf -> p_atm.

*/

#ifndef AMANZI_FLOWRELATIONS_ONE_SURF_KR_MODEL_
#define AMANZI_FLOWRELATIONS_ONE_SURF_KR_MODEL_

#include "Teuchos_ParameterList.hpp"
#include "surface_relperm_model.hh"
#include "dbc.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class OneSurfaceRelPermModel : public SurfaceRelPermModel {
 public:
  OneSurfaceRelPermModel(Teuchos::ParameterList& list) {}

  virtual bool TemperatureDependent() { return false; }

  virtual double SurfaceRelPerm(double uf, double h) { return 1.; }

  virtual double DSurfaceRelPermDUnfrozenFraction(double uf, double h) { return 0.; }

  virtual double DSurfaceRelPermDPondedDepth(double uf, double h) { return 0.; }

 private:
  static Utils::RegisteredFactory<SurfaceRelPermModel, OneSurfaceRelPermModel> reg_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
