/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Ugly hackjob to enable direct evaluation of the full model, on a single
  WRM/region.  This is bypassing much of the "niceness" of the framework, but
  seems necessary for solving a cell-wise correction equation.

  Uses intensive, not extensive, forms.

*/

#include "exceptions.hh"
#include "State.hh"

#include "eos_evaluator.hh"
#include "eos.hh"
#include "wrm_partition.hh"
#include "wrm_evaluator.hh"
#include "wrm.hh"
#include "molar_fraction_gas_evaluator.hh"
#include "vapor_pressure_relation.hh"
#include "pc_liquid_evaluator.hh"
#include "pc_liq_atm.hh"
#include "iem_evaluator.hh"
#include "iem.hh"
#include "iem_water_vapor_evaluator.hh"
#include "iem_water_vapor.hh"

#include "thermal_richards_model.hh"

namespace Amanzi {
namespace ATS_Physics {

#define DEBUG_FLAG 1

void
ThermalRichardsModel::InitializeModel(const Teuchos::Ptr<State>& S)
{
  // these are not yet initialized
  rho_rock_ = -1.;
  p_atm_ = -1.e12;

  // Grab the models.
  // get the WRM models and their regions
  Teuchos::RCP<Evaluator> me = S->GetEvaluator("saturation_gas");
  Teuchos::RCP<Flow::WRMEvaluator> wrm_me = Teuchos::rcp_dynamic_cast<Flow::WRMEvaluator>(me);
  AMANZI_ASSERT(wrm_me != Teuchos::null);
  Teuchos::RCP<Flow::WRMPartition> wrms = wrm_me->get_WRMs();

  // this needs fixed eventually, but for now assuming one WRM, and therefore
  // one model --etc
  AMANZI_ASSERT(wrms->second.size() == 1);

  // -- WRMs
  wrm_ = wrms->second[0];

  // -- liquid EOS
  me = S->GetEvaluator("molar_density_liquid");
  Teuchos::RCP<Relations::EOSEvaluator> eos_liquid_me =
    Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(me);
  AMANZI_ASSERT(eos_liquid_me != Teuchos::null);
  liquid_eos_ = eos_liquid_me->get_EOS();

  // -- gas EOS
  me = S->GetEvaluator("molar_density_gas");
  Teuchos::RCP<Relations::EOSEvaluator> eos_gas_me =
    Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(me);
  AMANZI_ASSERT(eos_gas_me != Teuchos::null);
  gas_eos_ = eos_gas_me->get_EOS();

  // -- gas vapor pressure
  me = S->GetEvaluator("mol_frac_gas");
  Teuchos::RCP<Relations::MolarFractionGasEvaluator> mol_frac_me =
    Teuchos::rcp_dynamic_cast<Relations::MolarFractionGasEvaluator>(me);
  AMANZI_ASSERT(mol_frac_me != Teuchos::null);
  vpr_ = mol_frac_me->get_VaporPressureRelation();

  // -- capillary pressure for liq/gas
  me = S->GetEvaluator("capillary_pressure_gas_liq");
  Teuchos::RCP<Flow::PCLiquidEvaluator> pc_liq_me =
    Teuchos::rcp_dynamic_cast<Flow::PCLiquidEvaluator>(me);
  AMANZI_ASSERT(pc_liq_me != Teuchos::null);
  pc_l_ = pc_liq_me->get_PCLiqAtm();

  // -- iem for liquid
  me = S->GetEvaluator("internal_energy_liquid");
  Teuchos::RCP<Energy::IEMEvaluator> iem_liquid_me =
    Teuchos::rcp_dynamic_cast<Energy::IEMEvaluator>(me);
  AMANZI_ASSERT(iem_liquid_me != Teuchos::null);
  liquid_iem_ = iem_liquid_me->get_IEM();

  // -- iem for gas
  me = S->GetEvaluator("internal_energy_gas");
  Teuchos::RCP<Energy::IEMWaterVaporEvaluator> iem_gas_me =
    Teuchos::rcp_dynamic_cast<Energy::IEMWaterVaporEvaluator>(me);
  AMANZI_ASSERT(iem_gas_me != Teuchos::null);
  gas_iem_ = iem_gas_me->get_IEM();

  // -- iem for rock
  me = S->GetEvaluator("internal_energy_rock");
  Teuchos::RCP<Energy::IEMEvaluator> iem_rock_me =
    Teuchos::rcp_dynamic_cast<Energy::IEMEvaluator>(me);
  AMANZI_ASSERT(iem_rock_me != Teuchos::null);
  rock_iem_ = iem_rock_me->get_IEM();
}


void
ThermalRichardsModel::UpdateModel(const Teuchos::Ptr<State>& S)
{
  // update scalars
  rho_rock_ = *S->GetScalarData("rho_rock");
  p_atm_ = *S->GetScalarData("atmospheric_pressure", Tags::DEFAULT);
  poro_ = (*S->Get<CompositeVector>("base_porosity").ViewComponent("cell"))[0][c];

  AMANZI_ASSERT(IsSetUp_());
}

// ----------------------------------------------------------------------
// Lightweight wrapper to forward-evaluate the model.
// ----------------------------------------------------------------------
int
ThermalRichardsModel::Evaluate(double T, double p, double poro, double& energy, double& wc)
{
  AmanziGeometry::Point res(2);
  int ierr = EvaluateEnergyAndWaterContent_(T, p, poro, res);
  energy = res[0];
  wc = res[1];
  return ierr;
};


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
ThermalRichardsModel::InverseEvaluate(double energy, double wc, double poro, double& T, double& p)
{
  // -- scaling for the norms
  double wc_scale = 10.;
  double e_scale = 10000.;
  double tol = 1.e-6;
  double max_steps = 100;
  double stepnum = 0;

  // get the initial residual
  AmanziGeometry::Point res(2);
  WhetStone::Tensor jac(2, 2);
  int ierr = EvaluateEnergyAndWaterContentAndJacobian_(T, p, poro, res, jac);
  if (ierr) {
    std::cout << "Error in evaluation: " << ierr << std::endl;
    return ierr + 10;
  }

#if DEBUG_FLAG
  std::cout << "Inverse Evaluating, e=" << energy << ", wc=" << wc << std::endl;
  std::cout << "   guess T,p (res) = " << T << ", " << p << " (" << res[0] << ", " << res[1] << ")"
            << std::endl;
#endif

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

    if (std::abs(detJ) < 1.e-12) {
      std::cout << " Zero determinate of Jacobian:" << std::endl;
      std::cout << "   [" << jac(0, 0) << "," << jac(0, 1) << "]" << std::endl;
      std::cout << "   [" << jac(1, 0) << "," << jac(1, 1) << "]" << std::endl;
      std::cout << "  at T,p = " << x_tmp[0] << ", " << x_tmp[1] << std::endl;
      std::cout << "  with res(e,wc) = " << res[0] << ", " << res[1] << std::endl;
      return 1;

    } else {
      jac.Inverse();
      correction = jac * res;
    }


    // perform the update
    x_tmp = x - correction;
    ierr = EvaluateEnergyAndWaterContentAndJacobian_(x_tmp[0], x_tmp[1], poro, res, jac);
    if (ierr) {
      std::cout << "Error in evaluation: " << ierr << std::endl;
      return ierr + 10;
    }
    res = res - f;

#if DEBUG_FLAG
    std::cout << "  Iter: " << stepnum;
    std::cout << " corrected T,p (res) = " << x_tmp[0] << ", " << x_tmp[1] << " (" << res[0] << ", "
              << res[1] << ")" << std::endl;
#endif

    // check convergence and damping
    scaled_res[0] = res[0] / e_scale;
    scaled_res[1] = res[1] / wc_scale;
    double norm_new = AmanziGeometry::norm(scaled_res);

    double damp = 1.;
    bool backtracking_required = false;
    while (norm_new > norm) {
      backtracking_required = true;

      // backtrack
      damp *= 0.5;
      x_tmp = x - (damp * correction);

      // evaluate the damped value
      ierr = EvaluateEnergyAndWaterContent_(x_tmp[0], x_tmp[1], poro, res);
      if (ierr) {
        std::cout << "Error in evaluation: " << ierr << std::endl;
        return ierr + 10;
      }
      res = res - f;

#if DEBUG_FLAG
      std::cout << "    Damping: " << stepnum;
      std::cout << " corrected T,p (res) = " << x_tmp[0] << ", " << x_tmp[1] << " (" << res[0]
                << ", " << res[1] << ")" << std::endl;
#endif

      // check the new residual
      scaled_res[0] = res[0] / e_scale;
      scaled_res[1] = res[1] / wc_scale;
      norm_new = AmanziGeometry::norm(scaled_res);
    }

    if (backtracking_required) {
      // must recalculate the Jacobian at the new value
      ierr = EvaluateEnergyAndWaterContentAndJacobian_(x_tmp[0], x_tmp[1], poro, res, jac);
      if (ierr) {
        std::cout << "Error in evaluation: " << ierr << std::endl;
        return ierr + 10;
      }
      res = res - f;
    }

    // iterate
    x = x_tmp;
    norm = norm_new;

    converged = norm < tol;
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


bool
ThermalRichardsModel::IsSetUp_()
{
  if (wrm_ == Teuchos::null) return false;
  if (liquid_eos_ == Teuchos::null) return false;
  if (gas_eos_ == Teuchos::null) return false;
  if (pc_l_ == Teuchos::null) return false;
  if (vpr_ == Teuchos::null) return false;
  if (liquid_iem_ == Teuchos::null) return false;
  if (gas_iem_ == Teuchos::null) return false;
  if (rock_iem_ == Teuchos::null) return false;
  if (rho_rock_ < 0.) return false;
  if (p_atm_ < -1.e10) return false;
  return true;
}


int
ThermalRichardsModel::EvaluateEnergyAndWaterContent_(double T,
                                                     double p,
                                                     double poro,
                                                     AmanziGeometry::Point& result)
{
  int ierr = 0;
  std::vector<double> eos_param(2);

  try {
    double eff_p = std::max(p_atm_, p);

    eos_param[0] = T;
    eos_param[1] = eff_p;

    double rho_l = liquid_eos_->MolarDensity(eos_params);
    double mass_rho_l = liquid_eos_->MassDensity(eos_params);
    double rho_g = gas_eos_->MolarDensity(eos_params);
    double omega = vpr_->SaturatedVaporPressure(T) / p_atm_;

    double pc_l = pc_l_->CapillaryPressure(p, p_atm_);

    double s_l = wrm_->saturation(pc_l);
    double s_g = 1.0 - s_l;

    double u_l = liquid_iem_->InternalEnergy(T);
    double u_g = gas_iem_->InternalEnergy(T, omega);

    double u_rock = rock_iem_->InternalEnergy(T);

    // water content
    result[1] = poro * (rho_l * s_l + rho_g * s_g * omega);

    // energy
    result[0] =
      poro * (u_l * rho_l * s_l + u_g * rho_g * s_g) + (1.0 - poro) * (rho_rock_ * u_rock);
  } catch (const Exceptions::Amanzi_exception& e) {
    if (e.what() == std::string("Cut timestep")) {
      ierr = 1;
    }
  }

  return ierr;
}


int
ThermalRichardsModel::EvaluateEnergyAndWaterContentAndJacobian_(double T,
                                                                double p,
                                                                double poro,
                                                                AmanziGeometry::Point& result,
                                                                WhetStone::Tensor& jac)
{
  return EvaluateEnergyAndWaterContentAndJacobian_FD_(T, p, poro, result, jac);
}


int
ThermalRichardsModel::EvaluateEnergyAndWaterContentAndJacobian_FD_(double T,
                                                                   double p,
                                                                   double poro,
                                                                   AmanziGeometry::Point& result,
                                                                   WhetStone::Tensor& jac)
{
  double eps_T = 1.e-7;
  double eps_p = 1.e-3;

  int ierr = EvaluateEnergyAndWaterContent_(T, p, poro, result);
  if (ierr) return ierr;

  AmanziGeometry::Point test(result);
  // d / dT
  ierr = EvaluateEnergyAndWaterContent_(T + eps_T, p, poro, test);
  if (ierr) return ierr;

  jac(0, 0) = (test[0] - result[0]) / (eps_T);
  jac(1, 0) = (test[1] - result[1]) / (eps_T);

  // d / dp
  ierr = EvaluateEnergyAndWaterContent_(T, p + eps_p, poro, test);
  if (ierr) return ierr;

  jac(0, 1) = (test[0] - result[0]) / (eps_p);
  jac(1, 1) = (test[1] - result[1]) / (eps_p);
  return 0;
}


} // namespace ATS_Physics
} // namespace Amanzi
