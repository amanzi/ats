/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
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
#include "wrm_permafrost_evaluator.hh"
#include "wrm_permafrost_model.hh"
#include "molar_fraction_gas_evaluator.hh"
#include "vapor_pressure_relation.hh"
#include "pc_liquid_evaluator.hh"
#include "pc_ice_evaluator.hh"
#include "pc_ice_water.hh"
#include "pc_liq_atm.hh"
#include "iem_evaluator.hh"
#include "iem.hh"
#include "compressible_porosity_evaluator.hh"
#include "compressible_porosity_model.hh"
#include "compressible_porosity_leijnse_evaluator.hh"
#include "compressible_porosity_leijnse_model.hh"
#include "liquid_ice_model.hh"

namespace Amanzi {

#define DEBUG_FLAG 0

void
LiquidIceModel::InitializeModel(const Teuchos::Ptr<State>& S,
                                const Tag& tag,
                                Teuchos::ParameterList& plist)
{
  tag_ = tag;

  // these are not yet initialized
  rho_rock_ = -1.;
  p_atm_ = -1.e12;
  domain = plist.get<std::string>("domain name", "");
  if (!domain.empty()) {
    mesh_ = S->GetMesh(domain);
  } else {
    mesh_ = S->GetMesh("domain");
  }

  Key liq_dens_key = Keys::readKey(plist, domain, "molar density liquid", "molar_density_liquid");
  Key ice_dens_key = Keys::readKey(plist, domain, "molar density ice", "molar_density_ice");
  Key iem_liq_key =
    Keys::readKey(plist, domain, "internal energy liquid", "internal_energy_liquid");
  Key iem_ice_key = Keys::readKey(plist, domain, "internal energy ice", "internal_energy_ice");
  Key iem_rock_key = Keys::readKey(plist, domain, "internal energy rock", "internal_energy_rock");
  Key si_key = Keys::readKey(plist, domain, "ice saturation", "saturation_ice");

  // Grab the models.
  // get the WRM models and their regions
  auto& sat_ice_eval = S->RequireEvaluator(si_key, tag);
  auto wrm_me = dynamic_cast<Flow::WRMPermafrostEvaluator*>(&sat_ice_eval);
  AMANZI_ASSERT(wrm_me != nullptr);
  wrms_ = wrm_me->get_WRMPermafrostModels();

  // -- liquid EOS
  auto& liq_eval = S->RequireEvaluator(liq_dens_key, tag);
  auto eos_liquid_me = dynamic_cast<Relations::EOSEvaluator*>(&liq_eval);
  AMANZI_ASSERT(eos_liquid_me != nullptr);
  liquid_eos_ = eos_liquid_me->get_EOS();

  // -- ice EOS
  auto& ice_eval = S->RequireEvaluator(ice_dens_key, tag);
  auto eos_ice_me = dynamic_cast<Relations::EOSEvaluator*>(&ice_eval);
  AMANZI_ASSERT(eos_ice_me != nullptr);
  ice_eos_ = eos_ice_me->get_EOS();

  // -- capillary pressure for ice/water
  use_pc_ice_ = plist.get<bool>("use pc_ice to determine sfc", true);
  if (use_pc_ice_) {
    auto& pc_ice_eval =
      S->RequireEvaluator(Keys::getKey(domain, "capillary_pressure_liq_ice"), tag);
    auto pc_ice_me = dynamic_cast<Flow::PCIceEvaluator*>(&pc_ice_eval);
    AMANZI_ASSERT(pc_ice_me != nullptr);
    pc_i_ = pc_ice_me->get_PCIceWater();
  }

  // -- capillary pressure for liq/gas
  auto& pc_gas_eval = S->RequireEvaluator(Keys::getKey(domain, "capillary_pressure_gas_liq"), tag);
  auto pc_liq_me = dynamic_cast<Flow::PCLiquidEvaluator*>(&pc_gas_eval);
  AMANZI_ASSERT(pc_liq_me != nullptr);
  pc_l_ = pc_liq_me->get_PCLiqAtm();

  // -- iem for liquid
  auto& iem_liq_eval = S->RequireEvaluator(iem_liq_key, tag);
  auto iem_liquid_me = dynamic_cast<Energy::IEMEvaluator*>(&iem_liq_eval);
  AMANZI_ASSERT(iem_liquid_me != nullptr);
  liquid_iem_ = iem_liquid_me->get_IEM();

  // -- iem for ice
  auto& iem_ice_eval = S->RequireEvaluator(iem_ice_key, tag);
  auto iem_ice_me = dynamic_cast<Energy::IEMEvaluator*>(&iem_ice_eval);
  AMANZI_ASSERT(iem_ice_me != nullptr);
  ice_iem_ = iem_ice_me->get_IEM();

  // -- iem for rock
  auto& iem_rock_eval = S->RequireEvaluator(iem_rock_key, tag);
  auto iem_rock_me = dynamic_cast<Energy::IEMEvaluator*>(&iem_rock_eval);
  AMANZI_ASSERT(iem_rock_me != nullptr);
  rock_iem_ = iem_rock_me->get_IEM();

  // -- porosity
  poro_leij_ = plist.get<bool>("porosity leijnse model", false);
  auto& poro_eval = S->RequireEvaluator(Keys::getKey(domain, "porosity"), tag);
  if (!poro_leij_) {
    auto poro_me = dynamic_cast<Flow::CompressiblePorosityEvaluator*>(&poro_eval);
    AMANZI_ASSERT(poro_me != nullptr);
    poro_models_ = poro_me->get_Models();
  } else {
    auto poro_me = dynamic_cast<Flow::CompressiblePorosityLeijnseEvaluator*>(&poro_eval);
    AMANZI_ASSERT(poro_me != nullptr);
    poro_leij_models_ = poro_me->get_Models();
  }
}


void
LiquidIceModel::UpdateModel(const Teuchos::Ptr<State>& S, int c)
{
  // update scalars
  p_atm_ = S->Get<double>("atmospheric_pressure", Tags::DEFAULT);
  rho_rock_ = (*S->Get<CompositeVector>(Keys::getKey(domain, "density_rock"), tag_)
                  .ViewComponent("cell"))[0][c];
  poro_ = (*S->Get<CompositeVector>(Keys::getKey(domain, "base_porosity"), tag_)
              .ViewComponent("cell"))[0][c];
  wrm_ = wrms_->second[(*wrms_->first)[c]];
  if (!poro_leij_)
    poro_model_ = poro_models_->second[(*poro_models_->first)[c]];
  else
    poro_leij_model_ = poro_leij_models_->second[(*poro_leij_models_->first)[c]];

  AMANZI_ASSERT(IsSetUp_());
}

bool
LiquidIceModel::IsSetUp_()
{
  if (wrm_ == Teuchos::null) return false;
  if (!poro_leij_) {
    if (poro_model_ == Teuchos::null) return false;
  } else {
    if (poro_leij_model_ == Teuchos::null) return false;
  }
  if (liquid_eos_ == Teuchos::null) return false;
  if (ice_eos_ == Teuchos::null) return false;
  if (pc_l_ == Teuchos::null) return false;
  if (liquid_iem_ == Teuchos::null) return false;
  if (ice_iem_ == Teuchos::null) return false;
  if (rock_iem_ == Teuchos::null) return false;
  if (rho_rock_ < 0.) return false;
  if (p_atm_ < -1.e10) return false;
  if (use_pc_ice_) {
    if (pc_i_ == Teuchos::null) return false;
  }
  return true;
}


bool
LiquidIceModel::Freezing(double T, double p)
{
  double eff_p = std::max(p_atm_, p);
  std::vector<double> eos_param(2);
  eos_param[0] = T;
  eos_param[1] = eff_p;

  double pc_l = pc_l_->CapillaryPressure(p, p_atm_);
  double pc_i;
  if (use_pc_ice_) {
    if (pc_i_->IsMolarBasis()) {
      double rho_l = liquid_eos_->MolarDensity(eos_param);
      pc_i = pc_i_->CapillaryPressure(T, rho_l);
    } else {
      double mass_rho_l = liquid_eos_->MassDensity(eos_param);
      pc_i = pc_i_->CapillaryPressure(T, mass_rho_l);
    }
  }

  return wrm_->freezing(T, pc_l, pc_i);
}


int
LiquidIceModel::EvaluateSaturations(double T, double p, double& s_gas, double& s_liq, double& s_ice)
{
  int ierr = 0;
  try {
    double eff_p = std::max(p_atm_, p);

    std::vector<double> eos_param(2);
    eos_param[0] = T;
    eos_param[1] = eff_p;

    double pc_l = pc_l_->CapillaryPressure(p, p_atm_);
    double pc_i;
    if (use_pc_ice_) {
      if (pc_i_->IsMolarBasis()) {
        double rho_l = liquid_eos_->MolarDensity(eos_param);
        pc_i = pc_i_->CapillaryPressure(T, rho_l);
      } else {
        double mass_rho_l = liquid_eos_->MassDensity(eos_param);
        pc_i = pc_i_->CapillaryPressure(T, mass_rho_l);
      }
    }

    double sats[3];
    if (use_pc_ice_) {
      wrm_->saturations(pc_l, pc_i, sats);
    } else {
      wrm_->saturations(pc_l, T, sats);
    }
    s_gas = sats[0];
    s_liq = sats[1];
    s_ice = sats[2];

  } catch (const Exceptions::Amanzi_exception& e) {
    if (e.what() == std::string("Cut timestep")) { ierr = 1; }
  }

  return ierr;
}

int
LiquidIceModel::EvaluateEnergyAndWaterContent_(double T, double p, AmanziGeometry::Point& result)
{
  if (T < 100.0 || T > 373.0) {
    return 1; // invalid temperature
  }
  int ierr = 0;
  try {
    double poro;
    if (!poro_leij_)
      poro = poro_model_->Porosity(poro_, p, p_atm_);
    else
      poro = poro_leij_model_->Porosity(poro_, p, p_atm_);

    double eff_p = std::max(p_atm_, p);

    std::vector<double> eos_param(2);
    eos_param[0] = T;
    eos_param[1] = eff_p;

    double rho_l = liquid_eos_->MolarDensity(eos_param);
    double rho_i = ice_eos_->MolarDensity(eos_param);

    double pc_i;
    if (use_pc_ice_) {
      if (pc_i_->IsMolarBasis()) {
        pc_i = pc_i_->CapillaryPressure(T, rho_l);
      } else {
        double mass_rho_l = liquid_eos_->MassDensity(eos_param);
        pc_i = pc_i_->CapillaryPressure(T, mass_rho_l);
      }
    }

    double pc_l = pc_l_->CapillaryPressure(p, p_atm_);

    double sats[3];
    if (use_pc_ice_) {
      wrm_->saturations(pc_l, pc_i, sats);
    } else {
      wrm_->saturations(pc_l, T, sats);
    }
    double s_g = sats[0];
    double s_l = sats[1];
    double s_i = sats[2];

    double u_l = liquid_iem_->InternalEnergy(T);
    double u_i = ice_iem_->InternalEnergy(T);

    double u_rock = rock_iem_->InternalEnergy(T);

    // water content
    result[1] = poro * (rho_l * s_l + rho_i * s_i);

    // energy
    result[0] =
      poro * (u_l * rho_l * s_l + u_i * rho_i * s_i) + (1.0 - poro_) * (rho_rock_ * u_rock);
  } catch (const Exceptions::Amanzi_exception& e) {
    if (e.what() == std::string("Cut timestep")) { ierr = 1; }
  }

  return ierr;
}

} // namespace Amanzi
