/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include "UnitTest++.h"

#include "wrm_van_genuchten.hh"
#include "wrm_old_permafrost_model.hh"
#include "pc_ice_water.hh"


TEST(permafrost)
{
  using namespace Amanzi::Flow::Flow;

  double m = 0.8;
  double alpha = 1.5e-4;
  double sr = 0.;
  double p_atm = 101325.;

  Teuchos::ParameterList plist;
  plist.set("van Genuchten m", m);
  plist.set("van Genuchten alpha", alpha);
  plist.set("van Genuchten residual saturation", sr);
  plist.set("van Genuchten smoothing interval width", 0.0);
  Teuchos::RCP<WRMVanGenuchten> wrm = Teuchos::rcp(new WRMVanGenuchten(plist));

  Teuchos::ParameterList plist3;
  PCIceWater pcice(plist3);
  double rho = 1000.;

  Teuchos::ParameterList plist2;
  WRMOldPermafrostModel p1(plist2);
  p1.set_WRM(wrm);

  double sats[3];
  double sats2[3];
  double sats3[3];

  // unsaturated, above freezing
  // -- value
  double pc_ice = pcice.CapillaryPressure(275., rho);
  double pc_liq = p_atm - 100000.;
  p1.saturations(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[1], wrm->saturation(pc_liq), 1.e-8);

  // -- derivatives
  double eps = 1.e-6;
  p1.saturations(pc_liq + eps, pc_ice, sats2);
  p1.saturations(pc_liq - eps, pc_ice, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2 * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2 * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2 * eps);
  p1.dsaturations_dpc_liq(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3);


  // saturated, above freezing
  // -- value
  pc_liq = -1000.;
  p1.saturations(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[1], wrm->saturation(pc_liq), 1.e-8);

  // -- derivatives
  p1.saturations(pc_liq + eps, pc_ice, sats2);
  p1.saturations(pc_liq - eps, pc_ice, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2 * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2 * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2 * eps);
  p1.dsaturations_dpc_liq(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3);


  // saturated, below freezing
  // -- value
  pc_ice = pcice.CapillaryPressure(273.1, rho);
  pc_liq = p_atm - 100000.;
  p1.saturations(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], 1.5071633730387701e-09, 1.e-12);
  CHECK_CLOSE(sats[1], 6.075025773954138e-06, 1.e-12);
  CHECK_CLOSE(sats[2], 0.9999939234670627, 1.e-8);

  // -- derivatives
  p1.saturations(pc_liq + eps, pc_ice, sats2);
  p1.saturations(pc_liq - eps, pc_ice, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2 * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2 * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2 * eps);
  p1.dsaturations_dpc_liq(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3 + 1.e-10);

  eps = 1.e-5;
  p1.saturations(pc_liq, pc_ice + eps, sats2);
  p1.saturations(pc_liq, pc_ice - eps, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2. * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2. * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2. * eps);
  p1.dsaturations_dpc_ice(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3 + 1.e-10);


  // unsaturated, below freezing
  // -- value
  pc_liq = p_atm - (-200000);
  p1.saturations(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], 0.962055704092613, 1.e-6);
  CHECK_CLOSE(sats[1], 2.305125759593384e-07, 1.e-12);
  CHECK_CLOSE(sats[2], 0.037944065394811, 1.e-10);

  // -- derivatives
  p1.saturations(pc_liq + eps, pc_ice, sats2);
  p1.saturations(pc_liq - eps, pc_ice, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2 * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2 * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2 * eps);
  p1.dsaturations_dpc_liq(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3);

  p1.saturations(pc_liq, pc_ice + eps, sats2);
  p1.saturations(pc_liq, pc_ice - eps, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2 * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2 * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2 * eps);
  p1.dsaturations_dpc_ice(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3);
}
