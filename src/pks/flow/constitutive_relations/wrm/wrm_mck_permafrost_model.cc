/*

McKenzie et al. (2007)'s soil freezing curve

 */

#include "dbc.hh"

#include "wrm.hh"
#include "wrm_mck_permafrost_model.hh"

namespace Amanzi {
namespace Flow {

// required methods from the base class
// sats[0] = sg, sats[1] = sl, sats[2] = si, S*(pc_liq) = sl + si
void
WRMMCKPermafrostModel::saturations(double pc_liq, double temperature,
        double (&sats)[3]) {
  double sr_ = wrm_->residualSaturation();
  AMANZI_ASSERT(sr_ >= 0.); AMANZI_ASSERT(sr_ < 1.);
  if (pc_liq <= 0.) { // saturated
    sats[0] = 0.; // liquid + ice == 1
    if (temperature >= T0_) {
      sats[2] = 0.;
      sats[1] = 1.;
    } else {
      sats[1] = sr_+(1.0-sr_)*std::exp(-std::pow((temperature-T0_)/w_,2));
      sats[2] = 1. - sats[1];
    }
  } else { // unsaturated
    sats[0] = 1. - wrm_->saturation(pc_liq); // liquid+ice=1-gas
    if (temperature >= T0_) {
      sats[2] = 0.;
      sats[1] = wrm_->saturation(pc_liq);
    } else {
      double s = wrm_->saturation(pc_liq);
      sats[1] = sr_+(s-sr_)*std::exp(-std::pow((temperature-T0_)/w_,2));
      sats[2] = s - sats[1];
    }
  }
}

void
WRMMCKPermafrostModel::dsaturations_dpc_liq(double pc_liq, double temperature,
        double (&dsats)[3]) {
  if (pc_liq <=0.) {
    dsats[0] = 0.;
    dsats[1] = 0.;
    dsats[2] = 0.;
  } else {
    dsats[0] = -wrm_->d_saturation(pc_liq);
    if (temperature >= T0_) {
      dsats[2] = 0.;
      dsats[1] = wrm_->d_saturation(pc_liq);
    } else {
      dsats[1] = wrm_->d_saturation(pc_liq)*std::exp(-std::pow((temperature-T0_)/w_,2));
      dsats[2] = wrm_->d_saturation(pc_liq)-dsats[1];
    }
  }
}


void
WRMMCKPermafrostModel::dsaturations_dpc_ice(double pc_liq, double temperature,
        double (&dsats)[3]) {
  double sr_ = wrm_->residualSaturation();
  AMANZI_ASSERT(sr_ >= 0.); AMANZI_ASSERT(sr_ < 1.);
  dsats[0] = 0.;
  if (temperature >= T0_) {
    dsats[2] = 0.;
    dsats[1] = 0.;
  } else {
    if (pc_liq <= 0.) {
      dsats[1] = 2*(1.-sr_)*std::exp(-std::pow((temperature-T0_)/w_,2))*(T0_-temperature)/std::pow(w_,2);
    } else {
      double s = wrm_->saturation(pc_liq);
      dsats[1] = 2*(s-sr_)*std::exp(-std::pow((temperature-T0_)/w_,2))*(T0_-temperature)/std::pow(w_,2);
    }
    dsats[2] = -dsats[1];
  }
}

} // namespace Flow
} // namespace Flow
