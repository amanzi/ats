/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

#include <cmath>
#include "dbc.hh"
#include "errors.hh"
#include "Spline.hh"

#include "wrm_brooks_corey.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMBrooksCorey::WRMBrooksCorey(Teuchos::ParameterList& plist) : plist_(plist)
{
  InitializeFromPlist_();
};


void
WRMBrooksCorey::InitializeFromPlist_()
{
  lambda_ = plist_.get<double>("Brooks-Corey lambda [-]");
  b_ = 1. / lambda_; // clapp hornberger b
  p_sat_ = plist_.get<double>("Brooks-Corey saturated matric suction [Pa]");
  sr_ = plist_.get<double>("residual saturation [-]", 0.0);

  s0_ = 1.0 - plist_.get<double>("smoothing interval width [saturation]", 0.0);
  if (s0_ < 1.) { fit_kr_.Setup(s0_, k_relative(s0_), d_k_relative(s0_), 1.0, 1.0, 0.0); }
}


/* ******************************************************************
 * Relative permeability formula: input is saturation.
 * The original curve is regulized on interval (s0, 1) using the
 * Hermite interpolant of order 3. Formulas (3.11)-(3.12).
 ****************************************************************** */
double
WRMBrooksCorey::k_relative(double s)
{
  if (s <= s0_) {
    double se = (s - sr_) / (1 - sr_);
    return pow(se, 2 * b_ + 3);
  } else if (s == 1.) {
    return 1.;
  } else {
    return fit_kr_(s);
  }
}


/* ******************************************************************
 * D Relative permeability / D saturation
 ****************************************************************** */
double
WRMBrooksCorey::d_k_relative(double s)
{
  if (s <= s0_) {
    double se = (s - sr_) / (1 - sr_);
    double dkdse = (2 * b_ + 3) * pow(se, 2 * b_ + 2);
    return dkdse / (1. - sr_);
  } else if (s == 1.) {
    return 0.0;
  } else {
    return fit_kr_.Derivative(s);
  }
}


/* ******************************************************************
 * Saturation formula (3.5)-(3.8).
 ****************************************************************** */
double
WRMBrooksCorey::saturation(double pc)
{
  if (pc <= p_sat_) {
    return 1.;
  } else {
    return pow(p_sat_ / pc, lambda_) * (1. - sr_) + sr_;
  }
}


/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 ****************************************************************** */
double
WRMBrooksCorey::d_saturation(double pc)
{
  if (pc <= p_sat_) {
    return 0.;
  } else {
    return -(1. - sr_) * lambda_ * pow(p_sat_ / pc, lambda_) / pc;
  }
}


/* ******************************************************************
 * Pressure as a function of saturation, formula (3.9).
 ****************************************************************** */
double
WRMBrooksCorey::capillaryPressure(double s)
{
  double se = (s - sr_) / (1. - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  return pow(se, -b_) * p_sat_;
}


/* ******************************************************************
 * Derivative of pressure formulat w.r.t. saturation.
 ****************************************************************** */
double
WRMBrooksCorey::d_capillaryPressure(double s)
{
  double se = (s - sr_) / (1. - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  return -b_ * p_sat_ / (1. - sr_) * pow(se, -b_ - 1.);
}

} // namespace Flow
} // namespace Amanzi
