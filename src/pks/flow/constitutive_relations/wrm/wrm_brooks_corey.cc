/*
  This is the flow component of the Amanzi code.

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <string>

#include "dbc.hh"
#include "wrm_brooks_corey.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

/* ******************************************************************
 * Setup fundamental parameters for this model.
 * Default value of the regularization interval is pc0 = 0.
 ****************************************************************** */
WRMBrooksCorey::WRMBrooksCorey(Teuchos::ParameterList& plist) :
    plist_(plist) {
  InitializeFromPlist_();
};

void WRMBrooksCorey::InitializeFromPlist_() {
  lambda_ = plist_.get<double>("Brooks Corey lambda");
  l_ = plist_.get<double>("Mualem exponent l", 0.5);
  alpha_ = plist_.get<double>("Brooks Corey alpha");
  sr_ = plist_.get<double>("residual saturation", 0.0);
  pc0_ = plist_.get<double>("smoothing interval width", 0.0);

  std::string fname = plist_.get<std::string>("Krel function name", "Mualem");
  if (fname == std::string("Mualem")) {
    function_ = FLOW_WRM_MUALEM;
    factor_ = -2.0 - (l_ + 2.0) * lambda_;
  } else if (fname == std::string("Burdine")) {
    function_ = FLOW_WRM_BURDINE;
    factor_ = -2.0 - 3.0 * lambda_;
  } else {
    ASSERT(0);
  }

  pc_bubble_ = 1.0 / alpha_;
  a_ = b_ = 0;

  if (pc0_ > 0.0) {
    double pc0_orig = pc0_;
    pc0_ += pc_bubble_;
    double k0 = k_relative(pc0_) - 1.0;
    double k0p = d_k_relative(pc0_);
    double pc0_2 = pc0_orig * pc0_orig;
    double pc0_3 = pc0_2 * pc0_orig;

    a_ = (3 * k0 - k0p * pc0_orig) / pc0_2;
    b_ = (k0p * pc0_orig - 2 * k0) / pc0_3;
  }
}


/* ******************************************************************
 * Relative permeability formula: input is capillary pressure pc.
 * The original curve is regulized on interval (0, pc0) using the
 * Hermite interpolant of order 3. Formulas (3.14)-(3.15).
 ****************************************************************** */
double WRMBrooksCorey::k_relative(double pc)
{
  if (pc <= pc_bubble_) {
    return 1.0;
  } else if (pc >= pc0_) {
    return pow(alpha_ * pc, factor_);
  } else {
    double pc_2 = (pc - pc_bubble_) * (pc - pc_bubble_);
    double pc_3 = pc_2 * (pc - pc_bubble_);
    return 1.0 + a_ * pc_2 + b_ * pc_3;
  }
}

/* ******************************************************************
 * Derivative of the original relative permeability w.r.t. capillary pressure.
 ****************************************************************** */
double WRMBrooksCorey::d_k_relative(double pc)
{
  if (pc > pc_bubble_) {
    return factor_ * alpha_ * pow(alpha_ * pc, factor_ - 1.0);
  } else {
    return 0.0;
  }
}


/* ******************************************************************
 * Saturation formula (3.5)-(3.8).
 ****************************************************************** */
double WRMBrooksCorey::saturation(double pc)
{
  if (pc > pc_bubble_) {
    return pow(alpha_ * pc, -lambda_) * (1.0 - sr_) + sr_;
  } else {
    return 1.0;
  }
}


/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 * Warning: remember that dSdP = -d_saturation.
 ****************************************************************** */
double WRMBrooksCorey::d_saturation(double pc)
{
  if (pc > pc_bubble_) {
    return -pow(alpha_ * pc, -lambda_ - 1.0) * (1.0 - sr_) * alpha_ * lambda_;
  } else {
    return 0.0;
  }
}


/* ******************************************************************
 * Pressure as a function of saturation, formula (3.9).
 ****************************************************************** */
double WRMBrooksCorey::capillaryPressure(double s)
{
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);

  return pow(se, -1.0/lambda_) / alpha_;
}


/* ******************************************************************
 * Pressure as a function of saturation, formula (3.9).
 ****************************************************************** */
double WRMBrooksCorey::d_capillaryPressure(double s)
{
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);

  return -1. / lambda_ * pow(se, -1.0/lambda_ - 1.) / alpha_ / (1. - sr_);
}


}  // namespace AmanziFlow
}
}
