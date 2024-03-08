/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*!
  Sakagucki-Zeng soil resistance model refered to Sakaguchi and Zeng (2009).
  Note that `"dessicated_zone_thickness`" is given by soil types.
  If it is not declared in `"WRM paramters`" through `"model parameters`"
  under state, default value 0.1 m is used for all soil types.
*/


#ifndef AMANZI_FLOWRELATIONS_SOIL_RESISTANCE_SAKAGUCKI_ZENG_MODEL_HH_
#define AMANZI_FLOWRELATIONS_SOIL_RESISTANCE_SAKAGUCKI_ZENG_MODEL_HH_

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"

namespace Amanzi {
namespace Flow {

class SoilResistanceSakaguckiZengModel {
 public:
  explicit SoilResistanceSakaguckiZengModel(Teuchos::ParameterList& plist) : plist_(plist)
  {
    InitializeFromPlist_();
  }


  double RsoilbySakagickiZeng(double sat_gas, double porosity)
  {
    double r_soil;
    if (sat_gas == 0.) {
      r_soil = 0.; // ponded water
    } else {
      double vp_diffusion = 2.2e-5 * std::pow(porosity, 2) * std::pow(1 - sr_, 2 + 3 * b_);
      double L_Rsoil = d_ * (std::exp(std::pow(sat_gas, 5)) - 1) / (std::exp(1) - 1);
      r_soil = L_Rsoil / vp_diffusion;
    }
    AMANZI_ASSERT(r_soil >= 0);
    return r_soil;
  }


  double DRsoilbySakagickiZengDSatGas(double sat_gas, double porosity)
  {
    double vp_diffusion = 2.2e-5 * std::pow(porosity, 2) * std::pow(1 - sr_, 2 + 3 * b_);
    double coef = d_ / (std::exp(1) - 1) / vp_diffusion;
    return coef * std::exp(std::pow(sat_gas, 5)) * 5 * std::pow(sat_gas, 4);
  }


  double DRsoilbySakagickiZengDPorosity(double sat_gas, double porosity)
  {
    double L_Rsoil = d_ * (std::exp(std::pow(sat_gas, 5)) - 1) / (std::exp(1) - 1);
    double coef = L_Rsoil / 2.2e-5 * std::pow(1 - sr_, -2 - 3 * b_);
    return -2 * coef * std::pow(porosity, -3);
  }

 protected:
  void InitializeFromPlist_()
  {
    d_ = plist_.get<double>("dessicated zone thickness [m]", 0.1);
    sr_ = plist_.get<double>("residual saturation [-]", 0.0);

    if (plist_.get<std::string>("wrm type") == "van Genuchten") {
      if (plist_.isParameter("van Genuchten m [-]")) {
        double m = plist_.get<double>("van Genuchten m [-]");
        double n = 1.0 / (1.0 - m);
        double lambda = (n - 1) * (1 - std::pow(0.5, n / (n - 1)));
        b_ = 1. / lambda;
      } else {
        double n = plist_.get<double>("van Genuchten n [-]");
        double lambda = (n - 1) * (1 - std::pow(0.5, n / (n - 1)));
        b_ = 1. / lambda;
      }
    } else if (plist_.get<std::string>("wrm type") == "Brooks-Corey") {
      double lambda = plist_.get<double>("Brooks-Corey lambda [-]");
      b_ = 1. / lambda;
    } else {
      b_ = plist_.get<double>("Clapp-Hornberger b [-]");
    }
  }

 protected:
  Teuchos::ParameterList plist_;
  double sr_;
  double d_;
  double b_;
};

} // namespace Flow
} // namespace Amanzi

#endif
