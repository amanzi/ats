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

#ifndef AMANZI_FLOWRELATIONS_ZERO_UNFROZEN_FRACTION_KR_MODEL_
#define AMANZI_FLOWRELATIONS_ZERO_UNFROZEN_FRACTION_KR_MODEL_

#include <cmath>
#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include "Factory.hh"
#include "surface_relperm_model.hh"

namespace Amanzi {
namespace Flow {

class ZeroUFRelPermModel : public SurfaceRelPermModel {
 public:
  ZeroUFRelPermModel(Teuchos::ParameterList& list);

  virtual bool TemperatureDependent() { return true; }

  virtual double SurfaceRelPerm(double uf, double h);

  virtual double DSurfaceRelPermDUnfrozenFraction(double uf, double h)
  {
    if(h<=0. || alpha_==0)return 0.;
    double factor=h>=h_cutoff_?1.:std::pow(std::sin(pi_*h/(2.*h_cutoff_)),2);
    return factor*alpha_*pi_/2.*std::cos(pi_*uf/2.)*std::pow(std::sin(pi_*uf/2.),alpha_-1);
  }

  virtual double DSurfaceRelPermDPondedDepth(double uf, double h)
  {
    if(h<=0. || h>=h_cutoff_)return 0.;
    return std::pow(std::sin(pi_*uf/2.),alpha_)*pi_/h_cutoff_*
      std::sin(pi_*h/(2.*h_cutoff_))*std::cos(pi_*h/(2.*h_cutoff_));
  }

 protected:
  Teuchos::ParameterList plist_;

  int alpha_; // must be an even integer
  const double pi_;
  double h_cutoff_;

 private:
  static Utils::RegisteredFactory<SurfaceRelPermModel, ZeroUFRelPermModel> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
