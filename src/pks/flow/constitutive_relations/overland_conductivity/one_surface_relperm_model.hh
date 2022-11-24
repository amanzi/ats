
/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  No special limits as p_surf -> p_atm.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_ONE_SURF_KR_MODEL_
#define AMANZI_FLOWRELATIONS_ONE_SURF_KR_MODEL_

#include "Teuchos_ParameterList.hpp"
#include "surface_relperm_model.hh"
#include "dbc.hh"
#include "Factory.hh"

namespace Amanzi {
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
} // namespace Amanzi

#endif
