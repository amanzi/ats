/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

Painter's permafrost model with freezing point depression.

 */

#include "dbc.hh"

#include "wrm.hh"
#include "wrm_fpd_permafrost_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// required methods from the base class
// sats[0] = sg, sats[1] = sl, sats[2] = si
void
WRMFPDPermafrostModel::saturations(double pc_liq, double pc_ice, double (&sats)[3])
{
  if (pc_liq <= 0.) {                   // saturated
    sats[0] = 0.;                       // gas
    sats[1] = wrm_->saturation(pc_ice); // liquid
    sats[2] = 1.0 - sats[1];            // ice
    AMANZI_ASSERT(sats[2] >= 0.);
  } else if (pc_ice <= pc_liq) {
    sats[2] = 0.;                       // ice
    sats[1] = wrm_->saturation(pc_liq); // liquid
    sats[0] = 1.0 - sats[1];            // gas
    AMANZI_ASSERT(sats[0] >= 0.);
  } else {
    sats[1] = wrm_->saturation(pc_ice);
    sats[2] = 1. - sats[1] / wrm_->saturation(pc_liq);
    sats[0] = 1. - sats[1] - sats[2];
    AMANZI_ASSERT(sats[2] >= 0.);
    AMANZI_ASSERT(sats[0] >= 0.);
  }
}

void
WRMFPDPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3])
{
  if (pc_liq <= 0.) { // saturated
    dsats[0] = 0.;    // gas
    dsats[1] = 0.;
    dsats[2] = 0.;
  } else if (pc_ice <= pc_liq) {
    dsats[2] = 0.;                         // ice
    dsats[1] = wrm_->d_saturation(pc_liq); // liquid
    dsats[0] = -dsats[1];                  // gas
  } else {
    dsats[1] = 0.;
    dsats[2] =
      wrm_->saturation(pc_ice) / std::pow(wrm_->saturation(pc_liq), 2) * wrm_->d_saturation(pc_liq);
    dsats[0] = -dsats[2];
  }
}


void
WRMFPDPermafrostModel::dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3])
{
  if (pc_liq <= 0.) {                      // saturated
    dsats[0] = 0.;                         // gas
    dsats[1] = wrm_->d_saturation(pc_ice); // liquid
    dsats[2] = -dsats[1];                  // ice
  } else if (pc_ice <= pc_liq) {
    dsats[2] = 0.; // ice
    dsats[1] = 0.;
    dsats[0] = 0.;
  } else {
    dsats[1] = wrm_->d_saturation(pc_ice);
    dsats[2] = -dsats[1] / wrm_->saturation(pc_liq);
    dsats[0] = -dsats[1] - dsats[2];
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
