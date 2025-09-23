/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

Sutra model for saturated

 */

#include "dbc.hh"

#include "wrm.hh"
#include "wrm_sutra_permafrost_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// required methods from the base class
// sats[0] = sg, sats[1] = sl, sats[2] = si
void
WRMSutraPermafrostModel::saturations(double pc_liq, double temperature, double (&sats)[3])
{
  // pc_ice is temperature
  double dT = T0_ - temperature;
  if (dT <= 0.) {
    sats[0] = 0.;
    sats[1] = 1.;
    sats[2] = 0.;
  } else if (dT > dT_) {
    sats[0] = 0.;
    sats[1] = sr_;
    sats[2] = 1. - sr_;
  } else {
    sats[0] = 0.;
    sats[1] = 1.0 - (1.0 - sr_) * (dT / dT_);
    sats[2] = 1.0 - sats[1];
  }
}

void
WRMSutraPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3])
{
  dsats[0] = 0.;
  dsats[1] = 0.;
  dsats[2] = 0.;
}


void
WRMSutraPermafrostModel::dsaturations_dpc_ice(double pc_liq, double temperature, double (&dsats)[3])
{
  double dT = T0_ - temperature;
  if (dT <= 0.) {
    dsats[0] = 0.;
    dsats[0] = 0.;
    dsats[1] = 0.;
  } else if (dT > dT_) {
    dsats[0] = 0.;
    dsats[1] = 0.;
    dsats[2] = 0.;
  } else {
    dsats[0] = 0.;
    dsats[1] = (1.0 - sr_) / dT_;
    dsats[2] = -dsats[1];
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
