/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

EWCModelBase provides some of the functionality of EWCModel for inverse
evaluating.

------------------------------------------------------------------------- */

#include "ewc_model_base.hh"

#define DEBUG_FLAG 0

namespace Amanzi {

// ----------------------------------------------------------------------
// Lightweight wrapper to forward-evaluate the model.
// ----------------------------------------------------------------------
int
EWCModelBase::Evaluate(double T, double p, double& energy, double& wc)
{
  AmanziGeometry::Point res(2);
  int ierr = EvaluateEnergyAndWaterContent_(T, p, res);
  energy = res[0];
  wc = res[1];
  return ierr;
};


/* ----------------------------------------------------------------------
Solves a given energy and water content (at a given, fixed porosity), for
temperature and pressure.

Note this cannot really work for a saturated cell, as d_wc / d{T,p} = 0
<
Error codes:

  1 = Singular Jacobian at some point in the evaluation.  Often this is
      because the cell is saturated or becomes saturated and below
      freezing.
  2 = Iteration did not converge in max_steps (hard-coded to be 100 for
      now).
---------------------------------------------------------------------- */
int
EWCModelBase::InverseEvaluate(double energy, double wc, double& T, double& p, bool verbose)
{
  // -- scaling for the norms
  double wc_scale = 1.;
  double e_scale = 1.;
  double T_corr_cap = 2.;
  double p_corr_cap = 200000.;
  double tol = 1.e-6;
  double max_steps = 100;
  double stepnum = 0;

  // get the initial residual
  AmanziGeometry::Point res(2);
  WhetStone::Tensor jac(2, 2);
  int ierr = EvaluateEnergyAndWaterContentAndJacobian_(T, p, res, jac);
  if (ierr) {
    std::cout << "Error in evaluation: " << ierr << std::endl;
    return ierr + 10;
  }

  if (verbose) {
    std::cout << "Inverse Evaluating, e=" << energy << ", wc=" << wc << std::endl;
    std::cout << "   guess T,p (res) = " << T << ", " << p << " (" << res[0] << ", " << res[1]
              << ")" << std::endl;
  }

  AmanziGeometry::Point f(2);
  f[0] = energy;
  f[1] = wc;
  res = res - f;

  // check convergence
  AmanziGeometry::Point scaled_res(res);
  scaled_res[0] = res[0] / e_scale;
  scaled_res[1] = res[1] / wc_scale;
  double norm = AmanziGeometry::norm(scaled_res);

  bool converged = norm < tol;

  // workspace
  AmanziGeometry::Point x(2);
  AmanziGeometry::Point x_tmp(2);
  x[0] = T;
  x[1] = p;
  x_tmp[0] = T;
  x_tmp[1] = p;

  while (!converged) {
    // calculate the update size
    double detJ = jac.Det();
    AmanziGeometry::Point correction;

    if (std::abs(detJ) < 1.e-20) {
      std::cout << " Zero determinant of Jacobian:" << std::endl;
      std::cout << "   [" << jac(0, 0) << "," << jac(0, 1) << "]" << std::endl;
      std::cout << "   [" << jac(1, 0) << "," << jac(1, 1) << "]" << std::endl;
      std::cout << "  at T,p = " << x_tmp[0] << ", " << x_tmp[1] << std::endl;
      std::cout << "  with res(e,wc) = " << res[0] << ", " << res[1] << std::endl;
      return 1;
    }


    jac.Inverse();
    correction = jac * res;

    // cap the correction
    double scale = 1.;
    if (std::abs(correction[0]) > T_corr_cap) {
      scale = T_corr_cap / std::abs(correction[0]);
    }
    if (std::abs(correction[1]) > p_corr_cap) {
      double pscale = p_corr_cap / std::abs(correction[1]);
      scale = std::min(scale, pscale);
    }
    correction *= scale;

    // perform the update
    x_tmp = x - correction;
    ierr = EvaluateEnergyAndWaterContentAndJacobian_(x_tmp[0], x_tmp[1], res, jac);
    if (ierr) {
      std::cout << "Error in evaluation: " << ierr << std::endl;
      return ierr + 10;
    }
    res = res - f;

    // check convergence and damping
    scaled_res[0] = res[0] / e_scale;
    scaled_res[1] = res[1] / wc_scale;
    double norm_new = AmanziGeometry::norm(scaled_res);

    if (verbose) {
      std::cout << "  Iter: " << stepnum;
      std::cout << " corrected T,p (res) [norm] = " << x_tmp[0] << ", " << x_tmp[1] << " ("
                << res[0] << ", " << res[1] << ") [" << norm_new << "]" << std::endl;
    }

    double damp = 1.;
    bool backtracking_required = false;
    while (norm_new > norm) {
      backtracking_required = true;

      // backtrack
      damp *= 0.5;
      x_tmp = x - (damp * correction);

      // evaluate the damped value
      ierr = EvaluateEnergyAndWaterContent_(x_tmp[0], x_tmp[1], res);
      if (ierr) {
        std::cout << "Error in evaluation: " << ierr << std::endl;
        return ierr + 10;
      }
      res = res - f;

      // check the new residual
      scaled_res[0] = res[0] / e_scale;
      scaled_res[1] = res[1] / wc_scale;
      norm_new = AmanziGeometry::norm(scaled_res);

      if (verbose) {
        std::cout << "    Damping: " << stepnum;
        std::cout << " corrected T,p (res) [norm] = " << x_tmp[0] << ", " << x_tmp[1] << " ("
                  << res[0] << ", " << res[1] << ") [" << norm_new << "]" << std::endl;
      }
    }

    if (backtracking_required) {
      // must recalculate the Jacobian at the new value
      ierr = EvaluateEnergyAndWaterContentAndJacobian_(x_tmp[0], x_tmp[1], res, jac);
      if (ierr) {
        std::cout << "Error in evaluation: " << ierr << std::endl;
        return ierr + 10;
      }
      res = res - f;
    }

    // iterate
    x = x_tmp;
    norm = norm_new;

    AmanziGeometry::Point scaled_correction = damp * correction;
    scaled_correction[1] = scaled_correction[1] / 100000.;
    converged = norm < tol || AmanziGeometry::norm(scaled_correction) < 1.e-10;

    stepnum++;
    if (stepnum > max_steps && !converged) {
      std::cout << " Nonconverged after " << max_steps << " steps with norm (tol) " << norm << " ("
                << tol << ")" << std::endl;
      return 2;
    }
  }

  T = x[0];
  p = x[1];
  return 0;
}


/* ----------------------------------------------------------------------
Solves a given energy and water content (at a given, fixed porosity), for
temperature and pressure.

Note this cannot really work for a saturated cell, as d_wc / d{T,p} = 0

Error codes:

  1 = Singular Jacobian at some point in the evaluation.  Often this is
      because the cell is saturated or becomes saturated and below
      freezing.
  2 = Iteration did not converge in max_steps (hard-coded to be 100 for
      now).
---------------------------------------------------------------------- */
int
EWCModelBase::InverseEvaluateEnergy(double energy, double p, double& T)
{
  // -- scaling for the norms
  double e_scale = 1.;
  double T_corr_cap = 2.;
  double tol = 1.e-6;
  double max_steps = 100;
  double stepnum = 0;

  // get the initial residual
  AmanziGeometry::Point res(2);
  WhetStone::Tensor jac(2, 2);
  int ierr = EvaluateEnergyAndWaterContentAndJacobian_(T, p, res, jac);
  if (ierr) {
    std::cout << "Error in evaluation: " << ierr << std::endl;
    return ierr + 10;
  }

#if DEBUG_FLAG
  std::cout << "Inverse Evaluating, e=" << energy << std::endl;
  std::cout << "   guess T,p (res) = " << T << ", " << p << " (" << res[0] << ")" << std::endl;
#endif


  // check convergence
  double f = res[0] - energy;
  double norm = std::abs(f);
  bool converged = norm < tol;

  // workspace
  double T_tmp = T;
  double T_tmp2 = T;

  while (!converged) {
    // calculate the update size
    double detJ = jac(0, 0);
    double correction;

    if (std::abs(detJ) < 1.e-20) {
      std::cout << " Zero determinant of Jacobian:" << std::endl;
      std::cout << "   [" << jac(0, 0) << "]" << std::endl;
      std::cout << "  at T,p = " << T_tmp2 << ", " << p << std::endl;
      std::cout << "  with res(e) = " << f << std::endl;
      return 1;
    }

    correction = f / detJ;

    // cap the correction
    if (std::abs(correction) > T_corr_cap) {
      correction = correction / std::abs(correction) * T_corr_cap;
    }

    // perform the update
    T_tmp2 = T_tmp - correction;
    ierr = EvaluateEnergyAndWaterContentAndJacobian_(T_tmp2, p, res, jac);
    if (ierr) {
      std::cout << "Error in evaluation: " << ierr << std::endl;
      return ierr + 10;
    }
    f = res[0] - energy;

    // check convergence and damping
    double norm_new = std::abs(f);

#if DEBUG_FLAG
    std::cout << "  Iter: " << stepnum;
    std::cout << " corrected T,p (res) [norm] = " << T_tmp2 << ", " << p << " (" << f << ")"
              << std::endl;
#endif

    double damp = 1.;
    bool backtracking_required = false;
    while (norm_new > norm) {
      backtracking_required = true;

      // backtrack
      damp *= 0.5;
      T_tmp2 = T_tmp - (damp * correction);

      // evaluate the damped value
      ierr = EvaluateEnergyAndWaterContent_(T_tmp2, p, res);
      if (ierr) {
        std::cout << "Error in evaluation: " << ierr << std::endl;
        return ierr + 10;
      }
      f = res[0] - energy;

      // check the new residual
      norm_new = std::abs(f);

#if DEBUG_FLAG
      std::cout << "    Damping: " << stepnum;
      std::cout << " corrected T,p (res) [norm] = " << T_tmp2 << ", " << p << " (" << f << ")"
                << std::endl;
#endif
    }

    if (backtracking_required) {
      // must recalculate the Jacobian at the new value
      ierr = EvaluateEnergyAndWaterContentAndJacobian_(T_tmp2, p, res, jac);
      if (ierr) {
        std::cout << "Error in evaluation: " << ierr << std::endl;
        return ierr + 10;
      }
      f = res[0] - energy;
    }

    // iterate
    T_tmp = T_tmp2;
    norm = norm_new;

    converged = norm < tol || std::abs(correction) < 1.e-3;
    stepnum++;
    if (stepnum > max_steps && !converged) {
      std::cout << " Nonconverged after " << max_steps << " steps with norm (tol) " << norm << " ("
                << tol << ")" << std::endl;
      return 2;
    }
  }

  T = T_tmp;
  return 0;
}


int
EWCModelBase::EvaluateEnergyAndWaterContentAndJacobian_(double T,
                                                        double p,
                                                        AmanziGeometry::Point& result,
                                                        WhetStone::Tensor& jac)
{
  return EvaluateEnergyAndWaterContentAndJacobian_FD_(T, p, result, jac);
}


int
EWCModelBase::EvaluateEnergyAndWaterContentAndJacobian_FD_(double T,
                                                           double p,
                                                           AmanziGeometry::Point& result,
                                                           WhetStone::Tensor& jac)
{
  double eps_T = 1.e-7;
  double eps_p = 1.e-3;

  int ierr = EvaluateEnergyAndWaterContent_(T, p, result);
  if (ierr) return ierr;
  AmanziGeometry::Point test(result);
  AmanziGeometry::Point test2(result);

  // d / dT
  jac(0, 0) = 0.;
  jac(1, 0) = 0.;

  bool done = false;
  int its = 0;
  while (!done) {
    ierr = EvaluateEnergyAndWaterContent_(T + eps_T, p, test);
    if (ierr) return ierr;

    jac(0, 0) = (test[0] - result[0]) / (eps_T);
    jac(1, 0) = (test[1] - result[1]) / (eps_T);

    its++;
    done = (std::abs(jac(0, 0)) > 1.e-12) || (std::abs(jac(1, 0)) > 1.e-12);
    done |= (its > 30);
    eps_T *= 2;
  }

  // d / dp
  jac(0, 1) = 0.;
  jac(1, 1) = 0.;

  // failure point seems to be d/dp = 0, and p seems to need to be centered
  done = false;
  its = 0;
  while (!done) {
    ierr = EvaluateEnergyAndWaterContent_(T, p + eps_p, test);
    if (ierr) return ierr;
    ierr = EvaluateEnergyAndWaterContent_(T, p - eps_p, test2);
    if (ierr) return ierr;

    jac(0, 1) = (test[0] - test2[0]) / (2 * eps_p);
    jac(1, 1) = (test[1] - test2[1]) / (2 * eps_p);

    its++;
    done = (std::abs(jac(0, 1)) > 1.e-12) || (std::abs(jac(1, 1)) > 1.e-12);
    done |= (its > 30);
    eps_p *= 2;
  }

  return 0;
}


} // namespace Amanzi
