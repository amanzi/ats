/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include "dbc.hh"
#include "errors.hh"
#include "Spline.hh"

#include "wrm_van_genuchten.hh"

namespace Amanzi {
namespace Flow {

const double FLOW_WRM_TOLERANCE = 1e-10;

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMVanGenuchten::WRMVanGenuchten(Teuchos::ParameterList& plist) : plist_(plist)
{
  InitializeFromPlist_();
};


/* ******************************************************************
* Relative permeability formula: input is liquid saturation.
* The original curve is regulized on interval (s0, 1) using the
* Hermite interpolant of order 3. Formulas (3.11)-(3.12).
****************************************************************** */
double
WRMVanGenuchten::k_relative(double s)
{
  if (s <= s0_) {
    double se = (s - sr_) / (1 - sr_);
    if (function_ == FLOW_WRM_MUALEM) {
      return pow(se, l_) * pow(1.0 - pow(1.0 - pow(se, 1.0 / m_), m_), 2.0);
    } else {
      return se * se * (1.0 - pow(1.0 - pow(se, 1.0 / m_), m_));
    }
  } else if (s == 1.0) {
    return 1.0;
  } else {
    return fit_kr_(s);
  }
}


/* ******************************************************************
 * D Relative permeability / D capillary pressure pc.
 ****************************************************************** */
double
WRMVanGenuchten::d_k_relative(double s)
{
  if (s <= s0_) {
    double se = (s - sr_) / (1 - sr_);

    double x = pow(se, 1.0 / m_);
    if (fabs(1.0 - x) < FLOW_WRM_TOLERANCE) return 0.0;
    if (fabs(x) < FLOW_WRM_TOLERANCE) return 0.0;

    double y = pow(1.0 - x, m_);
    double dkdse;
    if (function_ == FLOW_WRM_MUALEM)
      dkdse = (1.0 - y) * (l_ * (1.0 - y) + 2 * x * y / (1.0 - x)) * pow(se, l_ - 1.0);
    else
      dkdse = (2 * (1.0 - y) + x / (1.0 - x)) * se;

    return dkdse / (1 - sr_);

  } else if (s == 1.0) {
    return 0.0;
  } else {
    return fit_kr_.Derivative(s);
  }
}


/* ******************************************************************
 * Saturation formula (3.5)-(3.6).
 ****************************************************************** */
double
WRMVanGenuchten::saturation(double pc)
{
  if (pc > pc0_) {
    return std::pow(1.0 + std::pow(alpha_ * pc, n_), -m_) * (1.0 - sr_) + sr_;
  } else if (pc <= 0.) {
    return 1.0;
  } else {
    return fit_s_(pc);
  }
}


/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 ****************************************************************** */
double
WRMVanGenuchten::d_saturation(double pc)
{
  if (pc > pc0_) {
    return -m_ * n_ * std::pow(1.0 + std::pow(alpha_ * pc, n_), -m_ - 1.0) *
           std::pow(alpha_ * pc, n_ - 1) * alpha_ * (1.0 - sr_);
  } else if (pc <= 0.) {
    return 0.0;
  } else {
    return fit_s_.Derivative(pc);
  }
}

/* ******************************************************************
 * Pressure as a function of saturation.
 ****************************************************************** */
double
WRMVanGenuchten::capillaryPressure(double s)
{
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  if (se < 1.e-8) {
    return std::pow(se, -1.0 / (m_ * n_)) / alpha_;
  } else {
    return (std::pow(std::pow(se, -1.0 / m_) - 1.0, 1 / n_)) / alpha_;
  }
}


/* ******************************************************************
 * Derivative of pressure formulat w.r.t. saturation.
 ****************************************************************** */
double
WRMVanGenuchten::d_capillaryPressure(double s)
{
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  if (se < 1.e-8) {
    return -1.0 / (m_ * n_ * alpha_) * std::pow(se, -1.0 / (m_ * n_) - 1.) / (1.0 - sr_);
  } else {
    return -1.0 / (m_ * n_ * alpha_) * std::pow(std::pow(se, -1.0 / m_) - 1.0, 1 / n_ - 1.0) *
           std::pow(se, -1.0 / m_ - 1.0) / (1.0 - sr_);
  }
}


void
WRMVanGenuchten::InitializeFromPlist_()
{
  std::string fname = plist_.get<std::string>("Krel function name", "Mualem");
  if (fname == std::string("Mualem")) {
    function_ = FLOW_WRM_MUALEM;
  } else if (fname == std::string("Burdine")) {
    function_ = FLOW_WRM_BURDINE;
  } else {
    AMANZI_ASSERT(0);
  }

  // DEPRECATION ERROR
  if (plist_.isParameter("van Genuchten alpha")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"van Genuchten alpha\" is now \"van Genuchten alpha [Pa^-1]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist_.isParameter("van Genuchten n")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"van Genuchten n\" is now \"van Genuchten n [-]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist_.isParameter("van Genuchten m")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"van Genuchten m\" is now \"van Genuchten m [-]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist_.isParameter("residual saturation")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"residual saturation\" is now \"residual saturation [-]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist_.isParameter("Mualem exponent l")) {
    Errors::Message message(
      "WRM: DEPRECATION: parameter \"Mualem exponent l\" is now \"Mualem exponent l [-]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist_.isParameter("smoothing interval width")) {
    Errors::Message message("WRM: DEPRECATION: option \"smoothing interval width\" is now "
                            "\"smoothing interval width [saturation]\"");
    Exceptions::amanzi_throw(message);
  }
  if (plist_.isParameter("van Genuchten residual saturation")) {
    Errors::Message message("WRM: DEPRECATION: option \"van Genuchten residual saturation\" is now "
                            "\"residual saturation [-]\"");
    Exceptions::amanzi_throw(message);
  }

  alpha_ = plist_.get<double>("van Genuchten alpha [Pa^-1]");
  sr_ = plist_.get<double>("residual saturation [-]", 0.0);
  l_ = plist_.get<double>("Mualem exponent l [-]", 0.5);

  // map to n,m
  if (plist_.isParameter("van Genuchten m [-]")) {
    m_ = plist_.get<double>("van Genuchten m [-]");
    if (function_ == FLOW_WRM_MUALEM) {
      n_ = 1.0 / (1.0 - m_);
    } else {
      n_ = 2.0 / (1.0 - m_);
    }
  } else {
    n_ = plist_.get<double>("van Genuchten n [-]");
    if (function_ == FLOW_WRM_MUALEM) {
      m_ = 1.0 - 1.0 / n_;
    } else {
      m_ = 1.0 - 2.0 / n_;
    }
  }

  s0_ = 1.0 - plist_.get<double>("smoothing interval width [saturation]", 0.0);
  if (s0_ < 1.) { fit_kr_.Setup(s0_, k_relative(s0_), d_k_relative(s0_), 1.0, 1.0, 0.0); }

  pc0_ = plist_.get<double>("saturation smoothing interval [Pa]", 0.0);
  if (pc0_ > 0.) { fit_s_.Setup(0.0, 1.0, 0.0, pc0_, saturation(pc0_), d_saturation(pc0_)); }
};

/* ******************************************************************
* Suction formula: input is liquid saturation.
****************************************************************** */
double
WRMVanGenuchten::suction_head(double s)
{
  double se = (s - sr_) / (1 - sr_);
  if (se > FLOW_WRM_TOLERANCE) {
    return -(1. / alpha_) * pow(pow(se, -1. / m_) - 1, 1. - m_);
  } else {
    return -1e+10;
  }
}


/* ******************************************************************
 * D suction_head / D saturation
 ****************************************************************** */
double
WRMVanGenuchten::d_suction_head(double s)
{
  double se = (s - sr_) / (1 - sr_);

  double x = pow(se, -1.0 / m_);
  //  if (fabs(1.0 - x) < FLOW_WRM_TOLERANCE) return 0.0;

  double dpsidse;
  //if (function_ == FLOW_WRM_MUALEM)
  if (se > FLOW_WRM_TOLERANCE) {
    dpsidse = ((m_ - 1) / (alpha_ * m_ * se)) * pow(x - 1, -m_) * x;
  } else {
    dpsidse = 0.;
  }
  return -dpsidse / (1 - sr_);
}

} // namespace Flow
} // namespace Amanzi
