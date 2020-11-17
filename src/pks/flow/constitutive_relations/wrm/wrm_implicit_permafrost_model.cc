/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Painter's original, implicitly defined permafrost model.

#include <cmath>

#include "Epetra_SerialDenseMatrix.h"

#include "dbc.hh"
#include "errors.hh"
#include "Tensor.hh"
#include "Point.hh"

#include "wrm.hh"
#include "wrm_implicit_permafrost_model.hh"


namespace Amanzi {
namespace Flow {

// Constructor
WRMImplicitPermafrostModel::WRMImplicitPermafrostModel(Teuchos::ParameterList& plist) :
    WRMPermafrostModel(plist) {
  eps_ = plist_.get<double>("converged tolerance", 1.e-12);
  max_it_ = plist_.get<int>("max iterations", 100);
  deriv_regularization_ = plist_.get<double>("minimum dsi_dpressure magnitude", 1.e-10);
  solver_ = plist_.get<std::string>("solver algorithm [bisection/toms]", "bisection");
}

// Above freezing calculation methods:
// -- saturation calculation, above freezing
bool WRMImplicitPermafrostModel::sats_unfrozen_(double pc_liq,
        double pc_ice, double (&sats)[3]) {
  if (pc_ice <= 0.) {  // above freezing, s_i = 0, s_l/s_g by usual curve
    sats[2] = 0.; // ice
    sats[1] = wrm_->saturation(pc_liq);  // liquid
    sats[0] = 1.0 - sats[1];  // gas
    AMANZI_ASSERT(sats[0] >= 0.);
    return true;
  }
  return false;
}


// -- ds_dpcliq calculation, above freezing
bool WRMImplicitPermafrostModel::dsats_dpc_liq_unfrozen_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_ice <= 0.) {  // above freezing, s_i = 0, s_l/s_g by usual curve
    dsats[2] = 0.; // ice
    dsats[1] = wrm_->d_saturation(pc_liq);  // liquid
    dsats[0] = - dsats[1];  // gas
    return true;
  }
  return false;
}


// -- ds_dpcice calculation, above freezing
bool WRMImplicitPermafrostModel::dsats_dpc_ice_unfrozen_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_ice <= 0.) {  // above freezing, s_i = 0, s_l/s_g by usual curve
    dsats[2] = 0.; // ice
    dsats[1] = 0.;
    dsats[0] = 0.;
    return true;
  }
  return false;
}


// saturated calculation methods:
// -- saturation calculation, saturated
bool WRMImplicitPermafrostModel::sats_saturated_(double pc_liq,
        double pc_ice, double (&sats)[3]) {
  if (pc_liq <= 0.) { // fully saturated, s_g = 0, s_l/s_i by S_star(pc_ic)
    sats[0] = 0.; // gas
    sats[1] = wrm_->saturation(pc_ice); // liquid
    sats[2] = 1.0 - sats[1];  // ice
    AMANZI_ASSERT(sats[2] >= 0.);
    return true;
  }
  return false;
}


// -- ds_dpcliq calculation, saturated
bool WRMImplicitPermafrostModel::dsats_dpc_liq_saturated_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_liq <= 0.) { // fully saturated, s_g = 0, s_l/s_i by S_star(pc_ic)
    dsats[0] = 0.; // gas
    dsats[1] = 0.;
    dsats[2] = 0.;
    return true;
  }
  return false;
}


// -- ds_dpcice calculation, saturated
bool WRMImplicitPermafrostModel::dsats_dpc_ice_saturated_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  if (pc_liq <= 0.) { // fully saturated, s_g = 0, s_l/s_i by S_star(pc_ic)
    dsats[0] = 0.; // gas
    dsats[1] = wrm_->d_saturation(pc_ice); // liquid
    dsats[2] = - dsats[1];  // ice
    return true;
  }
  return false;
}


// partially frozen, unsaturated calculations
// -- saturation calculation, partially frozen, unsaturated
bool WRMImplicitPermafrostModel::sats_frozen_unsaturated_(double pc_liq,
        double pc_ice, double (&sats)[3]) {
  double si = si_frozen_unsaturated_(pc_liq, pc_ice);
  sats[2] = si;
  sats[1] = (1. - si) * wrm_->saturation(pc_liq);
  sats[0] = 1. - si - sats[1];
  AMANZI_ASSERT(sats[0] >= 0.);
  AMANZI_ASSERT(sats[1] >= 0.);
  AMANZI_ASSERT(sats[2] >= 0.);
  AMANZI_ASSERT(sats[0] <= 1.);
  AMANZI_ASSERT(sats[1] <= 1.);
  AMANZI_ASSERT(sats[2] <= 1.);
  return true;
}

// -- ds_dpcliq calculation, partially frozen, unsaturated
bool WRMImplicitPermafrostModel::dsats_dpc_liq_frozen_unsaturated_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  double si = si_frozen_unsaturated_(pc_liq, pc_ice);
  double dsi_dpcliq = dsi_dpc_liq_frozen_unsaturated_(pc_liq, pc_ice, si);
  dsats[2] = dsi_dpcliq;
  dsats[1] = (1. - si) * wrm_->d_saturation(pc_liq) - dsi_dpcliq * wrm_->saturation(pc_liq);
  dsats[0] = - dsats[1] - dsats[2];
  return true;
}

// -- ds_dpcice calculation, partially frozen, unsaturated
bool WRMImplicitPermafrostModel::dsats_dpc_ice_frozen_unsaturated_(double pc_liq,
        double pc_ice, double (&dsats)[3]) {
  double si = si_frozen_unsaturated_(pc_liq, pc_ice);
  double dsi_dpcice = dsi_dpc_ice_frozen_unsaturated_(pc_liq, pc_ice, si);
  dsats[2] = dsi_dpcice;
  dsats[1] = - dsi_dpcice * wrm_->saturation(pc_liq);
  dsats[0] = - dsats[1] - dsats[2];
  return true;
}


// -- si calculation, partially frozen, unsaturated
double WRMImplicitPermafrostModel::si_frozen_unsaturated_(double pc_liq, double pc_ice) {
  double si(0.);

  // check if we are in the splined region
  double cutoff(0.), si_cutoff(0.);
  DetermineSplineCutoff_(pc_liq, pc_ice, cutoff, si_cutoff);
  if (pc_liq > cutoff) {
    // outside of the spline
    si = si_frozen_unsaturated_nospline_(pc_liq, pc_ice);
  } else {
    // fit spline, evaluate
    double spline[4];
    FitSpline_(pc_ice, cutoff, si_cutoff, spline);
    si = ((spline[0] * pc_liq + spline[1]) * pc_liq + spline[2]) * pc_liq + spline[3];
    si = std::max(si, 0.);
    AMANZI_ASSERT(si <= 1.);
  }

  return si;
}


// -- dsi_dpcliq calculation, partially frozen, unsaturated
double WRMImplicitPermafrostModel::dsi_dpc_liq_frozen_unsaturated_(double pc_liq,
        double pc_ice, double si) {
  // check if we are in the splined region
  double cutoff(0.), si_cutoff(0.);
  double dsi(0.);

  DetermineSplineCutoff_(pc_liq, pc_ice, cutoff, si_cutoff);
  if (pc_liq > cutoff) {
    // outside of the spline
    dsi = dsi_dpc_liq_frozen_unsaturated_nospline_(pc_liq, pc_ice, si);
  } else {
    // fit spline, evaluate
    double spline[4];
    FitSpline_(pc_ice, cutoff, si_cutoff, spline);
    dsi = (3 * spline[0] * pc_liq + 2 * spline[1]) * pc_liq + spline[2];
  }

  // regularize
  if (std::abs(dsi) < deriv_regularization_) {
    dsi = dsi < 0. ? -deriv_regularization_ : dsi > 0. ? deriv_regularization_ : 0;
  }
  return dsi;
}


// -- dsi_dpcice calculation, partially frozen, unsaturated
double WRMImplicitPermafrostModel::dsi_dpc_ice_frozen_unsaturated_(double pc_liq,
        double pc_ice, double si) {
  // check if we are in the splined region
  double cutoff(0.), si_cutoff(0.);
  DetermineSplineCutoff_(pc_liq, pc_ice, cutoff, si_cutoff);
  if (pc_liq > cutoff) {
    // outside of the spline
    return dsi_dpc_ice_frozen_unsaturated_nospline_(pc_liq, pc_ice, si);
  } else {
    // fit two splines, difference neighboring splines
    double spline1[4], spline2[4];
    FitSpline_(pc_ice, cutoff, si_cutoff, spline1);

    double delta_pc_ice = std::max(.1, pc_ice / 100.);

    double pc_ice2 = pc_ice + delta_pc_ice;
    double si_cutoff2 = si_frozen_unsaturated_nospline_(cutoff, pc_ice2);
    FitSpline_(pc_ice2, cutoff, si_cutoff2, spline2);

    double dspline[4];
    dspline[0] = (spline2[0] - spline1[0]) / delta_pc_ice;
    dspline[1] = (spline2[1] - spline1[1]) / delta_pc_ice;
    dspline[2] = (spline2[2] - spline1[2]) / delta_pc_ice;
    dspline[3] = (spline2[3] - spline1[3]) / delta_pc_ice;

    return ((dspline[0] * pc_liq + dspline[1]) * pc_liq + dspline[2]) * pc_liq + dspline[3];
  }
}


// Helper methods for spline
// -- Determine the point beyond which the spline is not needed
bool WRMImplicitPermafrostModel::DetermineSplineCutoff_(double pc_liq, double pc_ice,
        double& cutoff, double& si) {
  cutoff = std::exp(std::floor(std::log(pc_liq)));
  bool done(false);
  while (!done) {
    try {
      si = si_frozen_unsaturated_nospline_(cutoff, pc_ice, true); // use the version that throws on error
    } catch (const Errors::CutTimeStep& e) {
      cutoff = std::exp(std::log(cutoff) + 1.);
      continue;
    }

    double sstar = wrm_->saturation(cutoff);
    if ((1. - si) * sstar + si < (1. - 1.e-16)) {
      done = true;
    } else {
      cutoff = std::exp(std::log(cutoff) + 1.);
    }
  }
  return true;
}


// -- Determine the coefficients of the spline
bool WRMImplicitPermafrostModel::FitSpline_(double pc_ice, double cutoff,
        double si_cutoff, double (&coefs)[4]) {
  double dsi_cutoff = dsi_dpc_liq_frozen_unsaturated_nospline_(cutoff, pc_ice, si_cutoff);


  // spline given by si = a * pc^3 + b * pc^3 + c * pc + d
  // constraint equations:
  //   1. si(0) = si_saturated  (right-continuous)
  //   2. si(cutoff) = si_cutoff (left-continuous)
  //   3. dsi_dpcliq(0) = dsi_saturated (continuous right deriv)
  //   4. dsi_dpcliq(cutoff) = dsi_cutoff (continuous left derivative)
  // solve for (a,b,c)

  double sats[3];
  sats_saturated_(-1.0, pc_ice, sats);
  double si_saturated = sats[2];
  double dsi_saturated = (si_cutoff - si_saturated) / cutoff;

  // form the linear system
  Amanzi::WhetStone::Tensor M(2,2);
  Amanzi::AmanziGeometry::Point rhs(2);

  // Equation 1:
  coefs[3] = si_saturated;

  // Equation 3
  coefs[2] = dsi_saturated;

  // Equation 2:
  M(0,0) = std::pow(cutoff,3);
  M(0,1) = std::pow(cutoff,2);
  rhs[0] = si_cutoff - (coefs[2] * cutoff + coefs[3]);

  // Equation 4
  M(1,0) = 3 * cutoff * cutoff;
  M(1,1) = 2 * cutoff;
  rhs[1] = dsi_cutoff - coefs[2];

  M.Inverse();
  coefs[0] = M(0,0) * rhs[0] + M(0,1) * rhs[1];
  coefs[1] = M(1,0) * rhs[0] + M(1,1) * rhs[1];
  return true;
}


// -- si calculation, outside of the splined region
double WRMImplicitPermafrostModel::si_frozen_unsaturated_nospline_(double pc_liq,
        double pc_ice, bool throw_ok) {
  // solve implicit equation for s_i
  SatIceFunctor_ func(pc_liq, pc_ice, wrm_);
  Tol_ tol(eps_);
  boost::uintmax_t max_it(max_it_);
  double left = 0.;
  double right = 1.;

  std::pair<double,double> result;
  try {
    if (solver_ == "bisection") {
      result = boost::math::tools::bisect(func, left, right, tol, max_it);
    } else if (solver_ == "toms") {
      result =
          boost::math::tools::toms748_solve(func, left, right, tol, max_it);
    } else {
      Errors::Message emsg("Unknown solver method");
      Exceptions::amanzi_throw(emsg);
    }
  } catch (const std::exception& e) {
    // this throw should not be caught
    std::stringstream estream;
    estream << "WRMImplicitPermafrostModel failed: " << e.what() << std::endl;
    Errors::Message emsg(estream.str());
    Exceptions::amanzi_throw(emsg);
  }

  double si = (result.first + result.second) / 2.;
  AMANZI_ASSERT(0. <= si && si <= 1.);

  if (max_it >= max_it_) {
    // did not converge?  May be ABS converged but not REL converged!
    if (!tol(func(si),0.)) {
      std::cout << "WRMImplicitPermafrostModel did not converge, " << max_it
                << " iterations, error = " << func(si) << ", s_i = " << si
                << ", PC_{lg,il} = " << pc_liq << "," << pc_ice << std::endl;
      if (throw_ok) {
        Exceptions::amanzi_throw(Errors::CutTimeStep());
      }
    }
  }
  return si;
}


// -- dsi_dpcliq calculation, outside of the splined region
double WRMImplicitPermafrostModel::dsi_dpc_liq_frozen_unsaturated_nospline_(double pc_liq,
        double pc_ice, double si) {
  // differentiate the implicit functor, solve for dsi_dpcliq
  double sstar =  wrm_->saturation(pc_liq);
  double sstarprime = wrm_->d_saturation(pc_liq);
  double tmp = (1.0 - si) * sstar;
  double tmpprime = (1.0 - si) * sstarprime;
  double G = - wrm_->d_saturation( pc_ice + wrm_->capillaryPressure( tmp + si))
      * wrm_->d_capillaryPressure( tmp + si );

  double numer = tmpprime * (1 + G);
  double denom = - sstar + G * (1.0 - sstar);

  return -numer / denom;
}


// -- dsi_pcice calculation, outside of the splined region
double WRMImplicitPermafrostModel::dsi_dpc_ice_frozen_unsaturated_nospline_(double pc_liq,
        double pc_ice, double si) {
  // differentiate the implicit functor, solve for dsi_dpcice
  double sstar =  wrm_->saturation(pc_liq);
  double tmp = (1.0 - si) * sstar;
  double G1 = wrm_->d_saturation( pc_ice + wrm_->capillaryPressure( tmp + si));
  double G2 = wrm_->d_capillaryPressure(tmp + si);

  return -G1 / (sstar + G1 * G2 * (1 - sstar));
}


// PUBLIC METHODS
// Calculate the saturation
void WRMImplicitPermafrostModel::saturations(double pc_liq, double pc_ice,
        double (&sats)[3]) {
  if (sats_unfrozen_(pc_liq, pc_ice, sats)) return;
  if (sats_saturated_(pc_liq, pc_ice, sats)) return;
  sats_frozen_unsaturated_(pc_liq, pc_ice, sats);
  return;
}

void WRMImplicitPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  if (dsats_dpc_liq_unfrozen_(pc_liq, pc_ice, dsats)) return;
  if (dsats_dpc_liq_saturated_(pc_liq, pc_ice, dsats)) return;
  dsats_dpc_liq_frozen_unsaturated_(pc_liq, pc_ice, dsats);
};

void WRMImplicitPermafrostModel::dsaturations_dpc_ice(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  if (dsats_dpc_ice_unfrozen_(pc_liq, pc_ice, dsats)) return;
  if (dsats_dpc_ice_saturated_(pc_liq, pc_ice, dsats)) return;
  dsats_dpc_ice_frozen_unsaturated_(pc_liq, pc_ice, dsats);
};


} // namespace
} // namespace
