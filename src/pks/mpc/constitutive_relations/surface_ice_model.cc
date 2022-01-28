/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Ugly hackjob to enable direct evaluation of the full model, on a single
  WRM/region.  This is bypassing much of the "niceness" of the framework, but
  seems necessary for solving a cell-wise correction equation.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "exceptions.hh"
#include "State.hh"

#include "eos_evaluator.hh"
#include "eos.hh"
#include "iem_evaluator.hh"
#include "iem.hh"
#include "unfrozen_fraction_evaluator.hh"
#include "unfrozen_fraction_model.hh"
#include "icy_height_evaluator.hh"
#include "icy_height_model.hh"

#include "surface_ice_model.hh"

namespace Amanzi {

void
SurfaceIceModel::InitializeModel(const Teuchos::Ptr<State>& S, Teuchos::ParameterList& plist)
{
  // NOTE: intentially using liquid instead of ice on the surface!
  Key domain = plist.get<std::string>("domain name", "surface");
  Key liq_dens_key = Keys::readKey(plist, domain, "molar density liquid", "molar_density_liquid");
  Key ice_dens_key = Keys::readKey(plist, domain, "molar density ice", "molar_density_ice"); // NOTE: NOT ICE!
  Key iem_liq_key = Keys::readKey(plist, domain, "internal energy liquid", "internal_energy_liquid");
  Key iem_ice_key = Keys::readKey(plist, domain, "internal energy ice", "internal_energy_ice");
  Key pd_key = Keys::readKey(plist, domain, "ponded depth", "ponded_depth");
  Key uf_key = Keys::readKey(plist, domain, "unfrozen fraction", "unfrozen_fraction");

  // these are not yet initialized
  M_ = 0.0180153;
  gz_ = -1.e12;
  p_atm_ = -1.e12;

  // Grab the models.
  // -- liquid EOS
  Teuchos::RCP<Evaluator> me = S->GetEvaluator(liq_dens_key);
  Teuchos::RCP<Relations::EOSEvaluator> eos_liquid_me =
      Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(me);
  AMANZI_ASSERT(eos_liquid_me != Teuchos::null);
  liquid_eos_ = eos_liquid_me->get_EOS();

  // -- ice EOS
  me = S->GetEvaluator(ice_dens_key);
  Teuchos::RCP<Relations::EOSEvaluator> eos_ice_me =
      Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(me);
  AMANZI_ASSERT(eos_ice_me != Teuchos::null);
  ice_eos_ = eos_ice_me->get_EOS();

  // -- iem for liquid
  me = S->GetEvaluator(iem_liq_key);
  Teuchos::RCP<Energy::IEMEvaluator> iem_liquid_me =
      Teuchos::rcp_dynamic_cast<Energy::IEMEvaluator>(me);
  AMANZI_ASSERT(iem_liquid_me != Teuchos::null);
  liquid_iem_ = iem_liquid_me->get_IEM();

  // -- iem for ice
  me = S->GetEvaluator(iem_ice_key);
  Teuchos::RCP<Energy::IEMEvaluator> iem_ice_me =
      Teuchos::rcp_dynamic_cast<Energy::IEMEvaluator>(me);
  AMANZI_ASSERT(iem_ice_me != Teuchos::null);
  ice_iem_ = iem_ice_me->get_IEM();

  // -- ponded depth evaluator
  me = S->GetEvaluator(pd_key);
  Teuchos::RCP<Flow::IcyHeightEvaluator> icy_h_me =
      Teuchos::rcp_dynamic_cast<Flow::IcyHeightEvaluator>(me);
  AMANZI_ASSERT(icy_h_me != Teuchos::null);
  pd_ = icy_h_me->get_IcyModel();

  // -- unfrozen fraction evaluator
  me = S->GetEvaluator(uf_key);
  Teuchos::RCP<Flow::UnfrozenFractionEvaluator> uf_me =
      Teuchos::rcp_dynamic_cast<Flow::UnfrozenFractionEvaluator>(me);
  AMANZI_ASSERT(uf_me != Teuchos::null);
  uf_ = uf_me->get_Model();
}

void
SurfaceIceModel::UpdateModel(const Teuchos::Ptr<State>& S, int c) {
  // update scalars
  p_atm_ = *S->GetScalarData("atmospheric_pressure", Tags::DEFAULT);
  gz_ = -((*S->GetConstantVectorData("gravity", Tags::DEFAULT))[2]);
  AMANZI_ASSERT(IsSetUp_());
}

bool
SurfaceIceModel::IsSetUp_() {
  if (pd_ == Teuchos::null) return false;
  if (uf_ == Teuchos::null) return false;
  if (liquid_eos_ == Teuchos::null) return false;
  if (ice_eos_ == Teuchos::null) return false;
  if (liquid_iem_ == Teuchos::null) return false;
  if (ice_iem_ == Teuchos::null) return false;
  if (p_atm_ < -1.e10) return false;
  if (gz_ < -1.e10) return false;
  if (M_ <= 0) return false;
  return true;
}

int
SurfaceIceModel::EvaluateEnergyAndWaterContent_(double T, double p,
        AmanziGeometry::Point& result) {
  if (T < 100) return 1; // invalid temperature
  int ierr = 0;
  std::vector<double> eos_param(2);

  try {
    // water content [mol / A]
    double WC = p < p_atm_ ? 0. : (p - p_atm_) / (gz_ * M_);

    // energy [J / A]
    // -- unfrozen fraction
    double uf = uf_->UnfrozenFraction(T);

    // -- densities
    eos_param[0] = T;
    eos_param[1] = p;
    double rho_l = liquid_eos_->MassDensity(eos_param);
    double n_l = rho_l / M_;
    double rho_i = ice_eos_->MassDensity(eos_param);
    double n_i = rho_i / M_;

    // -- ponded depth
    double h = pd_->Height(p, uf, rho_l, rho_i, p_atm_, gz_);

    // -- internal energies
    double u_l = liquid_iem_->InternalEnergy(T);
    double u_i = ice_iem_->InternalEnergy(T);

    // energy
    double E = h * ( uf * n_l * u_l + (1-uf) * n_i * u_i);

    // store solution
    result[1] = WC;
    result[0] = E;

  } catch (const Exceptions::Amanzi_exception& e) {
    if (e.what() == std::string("Cut time step")) {
      ierr = 1;
    }
  }

  return ierr;
}


} // namespace
