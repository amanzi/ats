/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*

Interfrost model for saturated

 */

#include "dbc.hh"

#include "wrm.hh"
#include "wrm_interfrost_permafrost_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// required methods from the base class
// sats[0] = sg, sats[1] = sl, sats[2] = si
void
WRMInterfrostPermafrostModel::saturations(double pc_liq, double temperature, double (&sats)[3])
{
  // pc_ice is temperature
  if (temperature >= 273.15) {
    sats[0] = 0.;
    sats[1] = 1.;
    sats[2] = 0.;
  } else {
    sats[0] = 0.;
    sats[1] = sr_ + (1 - sr_) * std::exp(-std::pow((temperature - 273.15) / W_, 2));
    sats[2] = 1.0 - sats[1];
  }
}

void
WRMInterfrostPermafrostModel::dsaturations_dpc_liq(double pc_liq, double pc_ice, double (&dsats)[3])
{
  dsats[0] = 0.;
  dsats[1] = 0.;
  dsats[2] = 0.;
}


void
WRMInterfrostPermafrostModel::dsaturations_dpc_ice(double pc_liq,
                                                   double temperature,
                                                   double (&dsats)[3])
{
  if (temperature >= 273.15) {
    dsats[0] = 0.;
    dsats[1] = 0.;
    dsats[2] = 0.;
  } else {
    dsats[0] = 0.;
    dsats[1] = (1 - sr_) * std::exp(-std::pow((temperature - 273.15) / W_, 2)) *
               (-2 * (temperature - 273.15) / W_) / W_;
    dsats[2] = -dsats[1];
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
