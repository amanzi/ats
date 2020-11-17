/*
Author: Ethan Coon

Painter's permafrost model with freezing point depression.

 */

#include "dbc.hh"

#include "wrm.hh"
#include "wrm_fpd_smoothed_permafrost_model.hh"

namespace Amanzi {
namespace Flow {

WRMFPDSmoothedPermafrostModel::WRMFPDSmoothedPermafrostModel(Teuchos::ParameterList& plist) :
    WRMPermafrostModel(plist)
{
  if (plist.isParameter("smoothing length [Pa]")) {
    dp_ = plist.get<double>("smoothing length [Pa]");
    AMANZI_ASSERT(dp_ > 0.);
  } else {
    double sigma_ice_liq = plist.get<double>("interfacial tension ice-water", 33.1);
    double sigma_gas_liq = plist.get<double>("interfacial tension air-water", 72.7);
    double T0 = plist.get<double>("reference temperature [K]", 273.15);
    double heat_fusion = plist.get<double>("latent heat [J kg^-1]", 3.34e5);
    double dens = plist.get<double>("water density [kg m^-3]", 999.87);
    double delT = plist.get<double>("smoothing width [K]", 1.0);
    dp_ = dens * sigma_gas_liq / sigma_ice_liq * heat_fusion * delT / T0;
    AMANZI_ASSERT(dp_ > 0.);
  }
}


// required methods from the base class
// sats[0] = sg, sats[1] = sl, sats[2] = si
void
WRMFPDSmoothedPermafrostModel::saturations(double pc_liq, double pc_ice,
        double (&sats)[3]) {
  if (pc_ice <= std::max(pc_liq,0.)) { // unfrozen
    sats[2] = 0.; // ice
    sats[1] = wrm_->saturation(pc_liq);  // liquid
    sats[0] = 1.0 - sats[1];  // gas

  } else {
    // freezing
    double sstar_ice = wrm_->saturation(pc_ice);
    double sstar_liq = wrm_->saturation(pc_liq);
    double sr = wrm_->residualSaturation();
    if (pc_liq <= 0.) { // saturated
      double sl_sm = (sstar_liq - sr)*std::exp( -pc_ice/dp_) + sr;
      sats[1] = std::max(sstar_ice, sl_sm);
      sats[2] = 1. - sats[1]; // ice
      sats[0] = 0.; // gas
    } else {
      double sl_sm = (sstar_liq - sr)*std::exp( (pc_liq - pc_ice) / dp_) + sr;
      sats[1] = std::max(sstar_ice, sl_sm);
      sats[2] = 1. - sats[1] / sstar_liq;
      sats[0] = 1. - sats[1] - sats[2];
    }
  }
}

void
WRMFPDSmoothedPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  if (pc_ice <= std::max(pc_liq,0.)) { // unfrozen
    dsats[2] = 0.; // ice
    dsats[1] = wrm_->d_saturation(pc_liq);  // liquid
    dsats[0] = - dsats[1];  // gas

  } else {
    // freezing
    double sstar_ice = wrm_->saturation(pc_ice);
    double sstar_liq = wrm_->saturation(pc_liq);
    double sr = wrm_->residualSaturation();

    if (pc_liq <= 0.) { // saturated
      double sl_sm = (sstar_liq - sr)*std::exp(-pc_ice/dp_) + sr;
      if (sstar_ice > sl_sm) { // max is given by sstar_ice
        dsats[1] = 0.;
        dsats[2] = 0.;
        dsats[0] = 0.;
      } else { // max is given by sl_sm
        double sstarprime_liq = wrm_->d_saturation(pc_liq);
        dsats[1] =  sstarprime_liq * std::exp(-pc_ice/dp_);
        dsats[2] = -dsats[1];
        dsats[0] = 0.; // no gas
      }

    } else { // unsaturated
      double sl_sm = (sstar_liq - sr)*std::exp( (pc_liq - pc_ice) / dp_) + sr;
      if (sstar_ice > sl_sm) { // max is given by sstar_ice
        dsats[1] = 0.;
        dsats[2] = sstar_ice / std::pow(sstar_liq,2)
            * wrm_->d_saturation(pc_liq);
        dsats[0] = - dsats[2];

      } else { // max is given by sl_sm
        double sstarprime_liq = wrm_->d_saturation(pc_liq);
        dsats[1] =  sstarprime_liq * std::exp( (pc_liq - pc_ice)/dp_ )
            + (sstar_liq - sr) * std::exp((pc_liq - pc_ice)/dp_) / dp_;
        dsats[2] = -dsats[1] / sstar_liq + sl_sm / std::pow(sstar_liq, 2) * sstarprime_liq;
        dsats[0] = -dsats[1] - dsats[2];
      }
    }
  }
}

void
WRMFPDSmoothedPermafrostModel::dsaturations_dpc_ice(double pc_liq, double pc_ice,
        double (&dsats)[3]) {
  if (pc_ice <= std::max(pc_liq,0.)) { // unfrozen
    dsats[0] = 0.;
    dsats[1] = 0.;
    dsats[2] = 0.;

  } else {
    // freezing
    double sstar_ice = wrm_->saturation(pc_ice);
    double sstar_liq = wrm_->saturation(pc_liq);
    double sr = wrm_->residualSaturation();

    if (pc_liq <= 0.) { // saturated
      double sl_sm = (sstar_liq - sr)*std::exp(-pc_ice/dp_) + sr;
      if (sstar_ice > sl_sm) { // max is given by sstar_ice
        dsats[1] = wrm_->d_saturation(pc_ice);
      } else { // max is given by sl_sm
        dsats[1] = -(sstar_liq - sr) * std::exp(-pc_ice/dp_) / dp_;
      }
      dsats[2] = -dsats[1];
      dsats[0] = 0.;

    } else { // unsaturated
      double sl_sm = (sstar_liq - sr)*std::exp( (pc_liq - pc_ice)/dp_) + sr;
      if (sstar_ice > sl_sm) { // max is given by sstar_ice
        dsats[1] = wrm_->d_saturation(pc_ice);
      } else { // max is given by sl_sm
        dsats[1] = -(sstar_liq - sr) * std::exp( (pc_liq - pc_ice)/dp_) / dp_;
      }
      dsats[2] = -dsats[1] / sstar_liq;
      dsats[0] = -dsats[1] - dsats[2];
    }
  }
}

} // namespace Flow
} // namespace Flow
