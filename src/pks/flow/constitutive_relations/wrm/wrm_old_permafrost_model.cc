/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

Painter's permafrost model.

 */

#include "wrm.hh"
#include "wrm_old_permafrost_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// sats[0] = s_g, sats[1] = s_l, sats[2] = s_i
void
WRMOldPermafrostModel::saturations(double pc_liq, double pc_ice, double (&sats)[3])
{
  if (pc_ice == 0.) {
    sats[2] = 0.;
    sats[1] = wrm_->saturation(pc_liq);
    sats[0] = 1.0 - sats[1];
  } else {
    double B = 1.0 / wrm_->saturation(pc_liq);
    double A = 1.0 / wrm_->saturation(pc_ice);

    sats[1] = 1.0 / (A + B - 1.);
    sats[0] = sats[1] * (B - 1.);
    sats[2] = sats[1] * (A - 1.);
  }
}


void
WRMOldPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3])
{
  if (pc_ice == 0.) {
    dsats[2] = 0.;
    dsats[1] = wrm_->d_saturation(pc_liq);
    dsats[0] = -dsats[1];
  } else {
    saturations(pc_liq, pc_ice, dsats);
    double sl = dsats[1];

    double B = 1.0 / wrm_->saturation(pc_liq);
    double A = 1.0 / wrm_->saturation(pc_ice);
    dsats[1] = sl * sl * B * B * wrm_->d_saturation(pc_liq);
    dsats[0] = dsats[1] * (B - 1.) + sl * (-B * B * wrm_->d_saturation(pc_liq));
    dsats[2] = dsats[1] * (A - 1.);
  }
}


void
WRMOldPermafrostModel::dsaturations_dpc_ice(double pc_liq, double pc_ice, double (&dsats)[3])
{
  if (pc_ice == 0.) {
    dsats[2] = 0.;
    dsats[1] = 0.;
    dsats[0] = 0.;
  } else {
    saturations(pc_liq, pc_ice, dsats);
    double sl = dsats[1];

    double B = 1.0 / wrm_->saturation(pc_liq);
    double A = 1.0 / wrm_->saturation(pc_ice);

    dsats[1] = sl * sl * A * A * wrm_->d_saturation(pc_ice);
    dsats[0] = dsats[1] * (B - 1.);
    dsats[2] = dsats[1] * (A - 1.) + sl * (-A * A * wrm_->d_saturation(pc_ice));
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
