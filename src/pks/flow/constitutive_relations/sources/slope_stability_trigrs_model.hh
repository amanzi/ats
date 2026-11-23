/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Bo Gao (gaob@ornl.gov)
*/

/*!
  TRIGRS slope stability model referred to Baum et al. (2008).
*/


#ifndef AMANZI_FLOW_SLOPE_STABILITY_TRIGRS_MODEL_HH_
#define AMANZI_FLOW_SLOPE_STABILITY_TRIGRS_MODEL_HH_

#include "Teuchos_ParameterList.hpp"
#include "dbc.hh"
#include <cmath>
#include <algorithm>

namespace Amanzi {
namespace Flow {

class SlopeStabilityTrigrsModel {
 public:
  explicit SlopeStabilityTrigrsModel(Teuchos::ParameterList& plist)
  {
    sr_ = plist.get<double>("residual saturation [-]", 0.0);
    AMANZI_ASSERT(sr_ >= 0. && sr_ < 1.);
    sigma_min_ = plist.get<double>("minimum slope [rad]", 1.e-6);
    stability_max_ = plist.get<double>("maximum slope stability [-]", 1.e10);
  }

  double SlopeStabilityByTrigrs(double phiprm, double sigma,
                                double cprm, double pressure,
                                double sat, double rhos,
                                double depth, double gravity,
                                double porosity, double rhow) const
  {
    // phiprm, sigma: rad; pressure: relative pressure
    if (sigma < sigma_min_) return stability_max_;
    double dry_ss = std::tan(phiprm) / std::tan(sigma);
    double numerator = cprm - Kappa_(sat) * pressure * std::tan(phiprm);
    double rho_bulk = BulkDensity_(rhos, porosity, rhow, sat);
    return dry_ss + numerator / (rho_bulk * Denom_(sigma, depth, gravity));
  }


  double DSlopeStabilityByTrigrsDpressure(double phiprm, double sigma,
                                          double sat, double rhos,
                                          double depth, double gravity,
                                          double porosity, double rhow) const
  {
    if (sigma < sigma_min_) return 0.;
    double rho_bulk = BulkDensity_(rhos, porosity, rhow, sat);
    return -Kappa_(sat) * std::tan(phiprm) / (rho_bulk * Denom_(sigma, depth, gravity));
  }


  double DSlopeStabilityByTrigrsDSat(double phiprm, double sigma,
                                     double cprm, double pressure, 
                                     double sat, double rhos,
                                     double depth, double gravity,
                                     double porosity, double rhow) const
  {
    if (sigma < sigma_min_) return 0.;

    double rho_bulk = BulkDensity_(rhos, porosity, rhow, sat);
    double coeff = -1. / std::pow(rho_bulk, 2.) / Denom_(sigma, depth, gravity);
    double rest0 = cprm * rhow * porosity;
    if (sat <= sr_) return coeff * rest0;

    double rest1 = BulkDensity_(rhos, porosity, rhow, sr_) * pressure * std::tan(phiprm) / (1. - sr_);
    return coeff * (rest0 + rest1);
  }


  double DSlopeStabilityByTrigrsDPorosity(double phiprm, double sigma,
                                          double cprm, double pressure, 
                                          double sat, double rhos,
                                          double depth, double gravity,
                                          double porosity, double rhow) const

  {
    if (sigma < sigma_min_) return 0.;

    double rho_bulk = BulkDensity_(rhos, porosity, rhow, sat);
    double coeff = (rhos - rhow * sat) / std::pow(rho_bulk, 2.) / Denom_(sigma, depth, gravity);
    double rest = cprm - Kappa_(sat) * pressure * std::tan(phiprm);
    return coeff * rest;
  }


 protected:
  double BulkDensity_(double rhos, double porosity, double rhow, double sat) const 
  {
    return rhos * (1. - porosity) + rhow * porosity * sat;
  }

  double Denom_(double sigma, double depth, double gravity) const
  {
    AMANZI_ASSERT(depth > 0.);
    return gravity * depth * std::sin(sigma) * std::cos(sigma);
  }

  double Kappa_(double sat) const
  {
    return std::clamp((sat - sr_) / (1.0 - sr_), 0., 1.);
  }

 protected:
  double sr_, sigma_min_, stability_max_;
};

} // namespace Flow
} // namespace Amanzi

#endif
