/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Scott Painter
           Daniil Svyatsky
*/

#include <cmath>
#include "dbc.hh"
#include "errors.hh"

#include "wrm_macropore.hh"

namespace Amanzi {
namespace Flow {

const double FLOW_WRM_TOLERANCE = 1e-10;

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMMacropore::WRMMacropore(Teuchos::ParameterList& plist)
  : plist_(plist)
{
  InitializeFromPlist_();
};


/* ******************************************************************
* Relative permeability formula: input is liquid saturation.
* The original curve is regulized on interval (s0, 1) using the
* Hermite interpolant of order 3. Formulas (3.11)-(3.12).
****************************************************************** */
double
WRMMacropore::k_relative(double s)
{
  return std::pow(s, a_);
}


/* ******************************************************************
 * D Relative permeability / D capillary pressure pc.
 ****************************************************************** */
double
WRMMacropore::d_k_relative(double s)
{
  if (std::abs(a_ - 1) < FLOW_WRM_TOLERANCE) {
    return 1.;
  } else {
    return a_ * std::pow(s, a_ - 1);
  }
}


/* ******************************************************************
 * Saturation formula
 ****************************************************************** */
double
WRMMacropore::saturation(double pc)
{
  if (pc <= 0) {
    return 1.0;
  } else if ((pc > 0) && (pc < Pwe_)) {
    return sr_ + (1 - sr_) * (Pwe_ - pc) / Pwe_;
  } else {
    return sr_ * std::exp((Pwe_ - pc) / dP_);
  }
}


/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 ****************************************************************** */
double
WRMMacropore::d_saturation(double pc)
{
  if (pc <= 0) {
    return 0.0;
  } else if ((pc > 0) && (pc < Pwe_)) {
    return -(1 - sr_) / Pwe_;
  } else {
    return -(sr_ / dP_) * std::exp((Pwe_ - pc) / dP_);
  }
}

/* ******************************************************************
 * Pressure as a function of saturation.
 ****************************************************************** */
double
WRMMacropore::capillaryPressure(double s)
{
  AMANZI_ASSERT(0);
  return 0.;
}


/* ******************************************************************
 * Derivative of pressure formulat w.r.t. saturation.
 ****************************************************************** */
double
WRMMacropore::d_capillaryPressure(double s)
{
  AMANZI_ASSERT(0);
  return 0.;
}


void
WRMMacropore::InitializeFromPlist_()
{
  a_ = plist_.get<double>("macropore exponent a [-]");
  sr_ = plist_.get<double>("residual saturation [-]");
  dP_ = plist_.get<double>("smoothing parameter [Pa]", 10);
  Pwe_ = plist_.get<double>("water entry pressure [Pa]");
};


} // namespace Flow
} // namespace Amanzi
