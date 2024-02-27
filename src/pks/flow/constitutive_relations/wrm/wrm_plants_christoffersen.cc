/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cmath>
#include <algorithm>

#include "dbc.hh"
#include "Brent.hh"

#include "wrm_plants_christoffersen.hh"

namespace Amanzi {
namespace Flow {

const double FLOW_WRM_TOLERANCE = 1e-10;

/* ******************************************************************
 * Setup fundamental parameters for this model.
 ****************************************************************** */
WRMPlantChristoffersen::WRMPlantChristoffersen(Teuchos::ParameterList& plist) : plist_(plist)
{
  InitializeFromPlist_();
};

double
WRMPlantChristoffersen::capillaryPressure(double s)
{
  return -potential(s) * 1E6;
}

double
WRMPlantChristoffersen::d_capillaryPressure(double s)
{
  return -d_potential(s) * 1E6;
}

double
WRMPlantChristoffersen::potential(double s) const
{
  double psi;
  if (s > sft_)
    psi = potentialLinear(s);
  else if (s > stlp_)
    psi = fmin(potentialSol(s) + potentialP(s), potentialLinear(s));
  else
    psi = potentialSol(s);
  return psi;
}

double
WRMPlantChristoffersen::d_potential(double s) const
{
  double dpsi;
  if (s > sft_)
    dpsi = d_potentialLinear(s);
  else if (s > stlp_) {
    if ((potentialSol(s) + potentialP(s)) < potentialLinear(s))
      dpsi = d_potentialSol(s) + d_potentialP(s);
    else
      dpsi = d_potentialLinear(s);
  } else
    dpsi = d_potentialSol(s);
  return dpsi;
}

double
WRMPlantChristoffersen::potentialLinear(double s) const
{
  return psi0_ - (mcap_ * (1 - s));
}


double
WRMPlantChristoffersen::d_potentialLinear(double s) const
{
  return mcap_;
}

double
WRMPlantChristoffersen::potentialSol(double s) const
{
  double sstar;
  sstar = (s - sr_) / (sft_ - sr_);
  return -std::abs(pi0_) / sstar;
}

double
WRMPlantChristoffersen::d_potentialSol(double s) const
{
  return (-std::abs(pi0_) * (sft_ - sr_) * sr_) / (std::pow((s - sr_), 2));
}

double
WRMPlantChristoffersen::potentialP(double s) const
{
  double sstar;
  sstar = (s - sr_) / (sft_ - sr_);
  return std::abs(pi0_) - (eps_ * (1 - sstar));
}

double
WRMPlantChristoffersen::d_potentialP(double s) const
{
  //double sstar;
  //sstar = (s - sr_) / (sft_ - sr_);
  //return std::abs(pi0_) - (eps_ * (1 - sstar));
  return eps_ / (sft_ - sr_);
}

double
WRMPlantChristoffersen::saturation(double pc)
{
  double psi;
  double se;

  psi = -pc * 1E-6;

  if (psi > psi0_)
    se = (1.0 - sr_) / (sft_ - sr_);
  else if (psi > psicapfttrans_)
    se = (1.0 - sr_) / (sft_ - sr_) - (psi0_ - psi) / mcapstar_;
  else if (psi > psitlp_) {
    double b = std::abs(pi0_) - psi - eps_;
    se = (-b + sqrt(std::pow(b, 2) + 4 * eps_ * std::abs(pi0_))) / (2 * eps_);
  } else
    se = -std::abs(pi0_) / psi;
  return (se * (sft_ - sr_)) + sr_;
}


double
WRMPlantChristoffersen::d_saturation(double pc)
{
  double psi;
  double dse;

  psi = -pc * 1E-6;
  double dpsi = -1E-6;

  if (psi > psi0_)
    dse = 0;
  else if (psi > psicapfttrans_)
    dse = dpsi / mcapstar_;
  else if (psi > psitlp_) {
    double b = std::abs(pi0_) - psi - eps_;
    double db = -dpsi;
    dse = (db + ((db * b) / sqrt(std::pow(b, 2) + 4 * eps_ * std::abs(pi0_)))) / (2 * eps_);
  } else
    dse = std::abs(pi0_) * (std::pow(psi, -2)) * (dpsi);
  return dse * (sft_ - sr_);
}

double
WRMPlantChristoffersen::k_relative(double pc)
{
  return (saturation(pc) - sr_) / (1 - sr_);
}


double
WRMPlantChristoffersen::d_k_relative(double pc)
{
  return (d_saturation(pc)) / (1 - sr_);
}

/* ******************************************************************
* Relative permeability formula: input is capillary pressure pc.
* The original curve is regulized on interval (0, pc0) using the
* Hermite interpolant of order 3. Formulas (3.11)-(3.12).
****************************************************************** */
/*
double WRMPlantChristoffersen::k_relative(double pc) {
  if (pc >= pc0_) {
    double se = std::pow(1.0 + pow(alpha_*pc, n_), -m_);
    if (function_ == FLOW_WRM_MUALEM) {
      return std::pow(se, l_) * pow(1.0 - pow(1.0 - pow(se, 1.0/m_), m_), 2.0);
    } else {
      return se * se * (1.0 - std::pow(1.0 - pow(se, 1.0/m_), m_));
    }
  } else if (pc <= 0.0) {
    return 1.0;
  } else {
    double pc_2 = pc * pc;
    double pc_3 = pc_2 * pc;
    return 1.0 + a_ * pc_2 + b_ * pc_3;
  }
}

*/
/* ******************************************************************
 * D Relative permeability / D capillary pressure pc.
 ****************************************************************** */
/*
double WRMPlantChristoffersen::d_k_relative(double pc) {
  if (pc >= pc0_) {
    double se = std::pow(1.0 + pow(alpha_*pc, n_), -m_);
    double dsdp = d_saturation(pc);

    double x = std::pow(se, 1.0 / m_);
    if (fstd::abs(1.0 - x) < FLOW_WRM_TOLERANCE) return 0.0;

    double y = std::pow(1.0 - x, m_);
    double dkdse;
    if (function_ == FLOW_WRM_MUALEM)
      dkdse = (1.0 - y) * (l_ * (1.0 - y) + 2 * x * y / (1.0 - x)) * std::pow(se, l_ - 1.0);
    else
      dkdse = (2 * (1.0 - y) + x / (1.0 - x)) * se;

    return dkdse * dsdp / (1 - sr_);

  } else if (pc <= 0.0) {
    return 0.0;
  } else {
    return 2*a_*pc + 3*b_*pc*pc;
  }
}

*/
/* ******************************************************************
 * Saturation formula (3.5)-(3.6).
 ****************************************************************** */
/*
double WRMPlantChristoffersen::saturation(double pc) {
  if (pc > 0.0) {
    return std::std::pow(1.0 + std::pow(alpha_*pc, n_), -m_) * (1.0 - sr_) + sr_;
  } else {
    return 1.0;
  }
}
*/

/* ******************************************************************
 * Derivative of the saturation formula w.r.t. capillary pressure.
 ****************************************************************** */
/*
double WRMPlantChristoffersen::d_saturation(double pc) {
  if (pc > 0.0) {
    return -m_*n_ * std::std::pow(1.0 + std::pow(alpha_*pc, n_), -m_-1.0) * std::pow(alpha_*pc, n_-1) * alpha_ * (1.0 - sr_);
  } else {
    return 0.0;
  }
}
*/
/* ******************************************************************
 * Pressure as a function of saturation.
 ****************************************************************** */
/*
double WRMPlantChristoffersen::capillaryPressure(double s) {
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  if (se < 1.e-8) {
    return std::std::pow(se, -1.0/(m_*n_)) / alpha_;
  } else {
    return (std::std::pow(std::pow(se, -1.0/m_) - 1.0, 1/n_)) / alpha_;
  }
}

*/
/* ******************************************************************
 * Derivative of pressure formulat w.r.t. saturation.
 ****************************************************************** */
/*
double WRMPlantChristoffersen::d_capillaryPressure(double s) {
  double se = (s - sr_) / (1.0 - sr_);
  se = std::min<double>(se, 1.0);
  se = std::max<double>(se, 1.e-40);
  if (se < 1.e-8) {
    return -1.0/(m_*n_*alpha_) * std::std::pow(se, -1.0/(m_*n_) - 1.)  / (1.0 - sr_);
  } else {
    return -1.0/(m_*n_*alpha_) * std::std::pow( std::pow(se, -1.0/m_) - 1.0, 1/n_-1.0)
        * std::std::pow(se, -1.0/m_ - 1.0) / (1.0 - sr_);
  }
}
*/

void
WRMPlantChristoffersen::InitializeFromPlist_()
{
  std::string fname = plist_.get<std::string>("Krel function name", "Mualem");
  if (fname == std::string("Mualem")) {
    function_ = FLOW_WRM_MUALEM;
  } else if (fname == std::string("Burdine")) {
    function_ = FLOW_WRM_BURDINE;
  } else {
    AMANZI_ASSERT(0);
  }

  /*
  alpha_ = plist_.get<double>("van Genuchten alpha");
  sr_ = plist_.get<double>("residual saturation", 0.0);
  l_ = plist_.get<double>("Mualem exponent l", 0.5);
  pc0_ = plist_.get<double>("smoothing interval width", 0.0);
  */

  sr_ = plist_.get<double>("residual saturation [-]", 0.0);
  stlp_ = plist_.get<double>("saturation at turgor loss [-]", 0.0);
  eps_ = plist_.get<double>("bulk elastic modulus [MPa]", 15.90459885634534);
  psi0_ = plist_.get<double>("water potential at full saturation [Pa]", -0.08);
  pi0_ = plist_.get<double>("osmotic potential at full turgor [Pa]", 0.0);
  psicap_ = plist_.get<double>("water potential at full turgor [Pa]", 0.0);
  double star = 1.0 - (std::abs(pi0_) / eps_);
  sft_ = (stlp_ - (1.0 - star) * sr_) / star;
  scap_ = 1.0 - 0.61 * (1.0 - stlp_);
  psitlp_ = -std::abs(pi0_) / ((stlp_ - sr_) / (sft_ - sr_));
  mcap_ = (psi0_ - psicap_) / (1.0 - scap_);
  mcapstar_ = mcap_ * (sft_ - sr_);

  /*
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

  a_ = b_ = 0.;
  if (pc0_ > 0) {
    double k0 = k_relative(pc0_) - 1.0;
    double k0p = d_k_relative(pc0_);
    double pc0_2 = pc0_ * pc0_;
    double pc0_3 = pc0_2 * pc0_;

    a_ = (3 * k0 - k0p * pc0_) / pc0_2;
    b_ = (k0p * pc0_ - 2 * k0) / pc0_3;

  }
  */

  int nits = 20;
  double tol(1.e-8);
  double result = Utils::findRootBrent(*this, 0.0, 1.0, tol, &nits);
  scapfttrans_ = result;
  psicapfttrans_ = potentialLinear(scapfttrans_);
  /*
  printf("sr: %f\nstlp: %f\neps: %f\npsi0: %f\npi0: %f\npsicap: %f\n", sr_, stlp_, eps_, psi0_, pi0_, psicap_);
  printf("sft: %f\nscap: %f\npsitlp: %f\nmcap: %f\nmcapstar: %f\nscapfttrans: %f\npsicapfttrans: %f\n",
          sft_, scap_, psitlp_, mcap_, mcapstar_, scapfttrans_, psicapfttrans_);
  printf("star: %f\n", star);
  */
};

} // namespace Flow
} // namespace Amanzi
