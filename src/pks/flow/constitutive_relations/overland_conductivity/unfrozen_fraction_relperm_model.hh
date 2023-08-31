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

#ifndef AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_KR_MODEL_
#define AMANZI_FLOWRELATIONS_UNFROZEN_FRACTION_KR_MODEL_

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "Factory.hh"
#include "surface_relperm_model.hh"

namespace Amanzi {
namespace Flow {

class UnfrozenFractionRelPermModel : public SurfaceRelPermModel {
 public:
  UnfrozenFractionRelPermModel(Teuchos::ParameterList& list);

  virtual bool TemperatureDependent() { return true; }

  virtual double SurfaceRelPerm(double uf, double h);

  virtual double DSurfaceRelPermDUnfrozenFraction(double uf, double h)
  {
    AMANZI_ASSERT(0);
    return 0.;
  }

  virtual double DSurfaceRelPermDPondedDepth(double uf, double h)
  {
    AMANZI_ASSERT(0);
    return 0.;
  }

 protected:
  Teuchos::ParameterList plist_;

  int alpha_; // must be an even integer
  const double pi_;

 private:
  static Utils::RegisteredFactory<SurfaceRelPermModel, UnfrozenFractionRelPermModel> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
