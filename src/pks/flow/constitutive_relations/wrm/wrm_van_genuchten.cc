/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov) 
  Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include "dbc.hh"
#include "Spline.hh"

#include "wrm_van_genuchten.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

const double FLOW_WRM_TOLERANCE = 1e-10;

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMVanGenuchten::WRMVanGenuchten(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};


/* ******************************************************************
* Relative permeability formula: input is liquid saturation.
* The original curve is regulized on interval (s0, 1) using the
* Hermite interpolant of order 3. Formulas (3.11)-(3.12).
****************************************************************** */
double WRMVanGenuchten::k_relative(double s) {
  if (s <= s0_) {
    double se = (s - sr_)/(1-sr_);
    if (function_ == FLOW_WRM_MUALEM) {
      return pow(se, l_) * pow(1.0 - pow(1.0 - pow(se, 1.0/m_), m_), 2.0);
    } else {
      return se * se * (1.0 - pow(1.0 - pow(se, 1.0/m_), m_));
    }
  } else if (s == 1.0) {
    return 1.0;
  } else {
    return fit_(s);
  }
}


/* ******************************************************************
 * D Relative permeability / D capillary pressure pc.
 ****************************************************************** */
double WRMVanGenuchten::d_k_relative(double s) {
  if (s <= s0_) {
    double se = (s - sr_)/(1-sr_);

    double x = pow(se, 1.0 / m_);
    if (fabs(1.0 - x) < FLOW_WRM_TOLERANCE) return 0.0;

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
    return fit_.Derivative(s);
  }
}


/* ******************************************************************
 * Saturation formula (3.5)-(3.6).
 ****************************************************************** */
double WRMVanGenuchten::saturation(double pc) {
  if (pc > 0.0) {
    return std::pow(1.0 + std::pow(alpha_*pc, n_), -m_) * (1.0 - sr_) + sr_;
  } else {
    return 1.0;
  }
}


/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 ****************************************************************** */
double WRMVanGenuchten::d_saturation(double pc) {
  if (pc > 0.0) {
    return -m_*n_ * std::pow(1.0 + std::pow(alpha_*pc, n_), -m_-1.0) * std::pow(alpha_*pc, n_-1) * alpha_ * (1.0 - sr_);
  } else {
    return 0.0;
  }
}

/* ******************************************************************
 * Pressure as a function of saturation.
 ****************************************************************** */
double WRMVanGenuchten::capillaryPressure(double s) {
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  if (se < 1.e-8) {
    return std::pow(se, -1.0/(m_*n_)) / alpha_;
  } else {
    return (std::pow(std::pow(se, -1.0/m_) - 1.0, 1/n_)) / alpha_;
  }
}


/* ******************************************************************
 * Derivative of pressure formulat w.r.t. saturation.
 ****************************************************************** */
double WRMVanGenuchten::d_capillaryPressure(double s) {
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  if (se < 1.e-8) {
    return -1.0/(m_*n_*alpha_) * std::pow(se, -1.0/(m_*n_) - 1.)  / (1.0 - sr_);
  } else {
    return -1.0/(m_*n_*alpha_) * std::pow( std::pow(se, -1.0/m_) - 1.0, 1/n_-1.0)
        * std::pow(se, -1.0/m_ - 1.0) / (1.0 - sr_);
  }
}


void WRMVanGenuchten::InitializeFromPlist_() {
  std::string fname = plist_.get<std::string>("Krel function name", "Mualem");
  if (fname == std::string("Mualem")) {
    function_ = FLOW_WRM_MUALEM;
  } else if (fname == std::string("Burdine")) {
    function_ = FLOW_WRM_BURDINE;
  } else {
    ASSERT(0);
  }

  alpha_ = plist_.get<double>("van Genuchten alpha");
  sr_ = plist_.get<double>("residual saturation", 0.0);
  l_ = plist_.get<double>("Mualem exponent l", 0.5);

  // map to n,m
  if (plist_.isParameter("van Genuchten m")) {
    m_ = plist_.get<double>("van Genuchten m");
    if (function_ == FLOW_WRM_MUALEM) {
      n_ = 1.0 / (1.0 - m_);
    } else {
      n_ = 2.0 / (1.0 - m_);
    }
  } else {
    n_ = plist_.get<double>("van Genuchten n");
    if (function_ == FLOW_WRM_MUALEM) {
      m_ = 1.0 - 1.0/n_;
    } else {
      m_ = 1.0 - 2.0/n_;
    }
  }

  s0_ = 1.0 - plist_.get<double>("smoothing interval width [saturation]", 0.0);
  if (s0_ < 1.) {
    fit_.Setup(s0_, k_relative(s0_), d_k_relative(s0_),
	       1.0, 1.0, 0.0);
  }  
  
};

}  // namespace
}  // namespace
}  // namespace
