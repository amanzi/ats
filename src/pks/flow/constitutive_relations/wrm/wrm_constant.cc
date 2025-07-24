/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <cmath>
#include "dbc.hh"
#include "errors.hh"
#include "Spline.hh"

#include "wrm_constant.hh"

namespace Amanzi {
namespace Flow {

const double FLOW_WRM_TOLERANCE = 1e-10;

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMConstant::WRMConstant(Teuchos::ParameterList& plist)
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
WRMConstant::k_relative(double s)
{
  return a_;
}


/* ******************************************************************
 * D Relative permeability / D capillary pressure pc.
 ****************************************************************** */
double
WRMConstant::d_k_relative(double s)
{
  return 0.;
}


/* ******************************************************************
 * Saturation formula
 ****************************************************************** */
double
WRMConstant::saturation(double pc)
{
  return sr_;
}


/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 ****************************************************************** */
double
WRMConstant::d_saturation(double pc)
{
  return 0.;
}

/* ******************************************************************
 * Pressure as a function of saturation.
 ****************************************************************** */
double
WRMConstant::capillaryPressure(double s)
{
  AMANZI_ASSERT(0);
  return 0;
}


/* ******************************************************************
 * Derivative of pressure formulat w.r.t. saturation.
 ****************************************************************** */
double
WRMConstant::d_capillaryPressure(double s)
{
  AMANZI_ASSERT(0);
  return 0;
}


void
WRMConstant::InitializeFromPlist_()
{
  a_ = plist_.get<double>("constant value [-]", 1.0);
  sr_ = plist_.get<double>("residual saturation [-]", 1.0);
};


} // namespace Flow
} // namespace Amanzi
