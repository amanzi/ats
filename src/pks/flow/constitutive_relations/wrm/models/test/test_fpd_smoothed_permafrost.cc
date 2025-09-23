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
#include "wrm_fpd_smoothed_permafrost_model.hh"
#include "pc_ice_water.hh"


TEST(fpdSmoothedPermafrost_Derivs)
{
  using namespace Amanzi::ATS_Physics::Flow;
  double eps = 1.e-4;

  double m = 0.8;
  double alpha = 1.5e-4;
  double sr = 0.1;
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
  WRMFPDSmoothedPermafrostModel p1(plist2);
  p1.set_WRM(wrm);

  double sats[3];
  double sats2[3];
  double sats3[3];

  // ensure the saturated, above freezing case reduces to the expected vals
  // -- value
  double pc_liq = -1000.;
  double pc_ice = pcice.CapillaryPressure(275., rho);

  p1.saturations(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], 0.0, 1.e-8);
  CHECK_CLOSE(sats[1], 1.0, 1.e-8);
  CHECK_CLOSE(sats[2], 0.0, 1.e-8);

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

  // ensure the unsaturated, above freezing case reduces to the WRM
  // -- value
  pc_ice = pcice.CapillaryPressure(275., rho);
  pc_liq = p_atm - 100000.;
  p1.saturations(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], 1. - wrm->saturation(pc_liq), 1.e-8);
  CHECK_CLOSE(sats[1], wrm->saturation(pc_liq), 1.e-8);
  CHECK_CLOSE(sats[2], 0., 1.e-8);

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

  // saturated, barely below freezing
  // -- value is not checkable without python implementation
  pc_ice = pcice.CapillaryPressure(273.1, rho);
  pc_liq = -5.;
  p1.saturations(pc_liq, pc_ice, sats);
  std::cout << "PC (i/l) = " << pc_ice << ", " << pc_liq << std::endl;
  std::cout << "  sats = " << sats[0] << ", " << sats[1] << ", " << sats[2] << std::endl;

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
  std::cout << "  dsat_dliq (FD) = " << sats2[0] << ", " << sats2[1] << ", " << sats2[2]
            << std::endl;
  std::cout << "  dsat_dliq = " << sats[0] << ", " << sats[1] << ", " << sats[2] << std::endl;

  p1.saturations(pc_liq, pc_ice + eps, sats2);
  p1.saturations(pc_liq, pc_ice - eps, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2. * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2. * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2. * eps);
  p1.dsaturations_dpc_ice(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3 + 1.e-10);
  std::cout << "  dsat_dice (FD) = " << sats2[0] << ", " << sats2[1] << ", " << sats2[2]
            << std::endl;
  std::cout << "  dsat_dice = " << sats[0] << ", " << sats[1] << ", " << sats[2] << std::endl;

  // saturated, well below freezing
  // -- value is not checkable without python implementation
  pc_ice = pcice.CapillaryPressure(265., rho);
  pc_liq = -500.;
  p1.saturations(pc_liq, pc_ice, sats);
  std::cout << "PC (i/l) = " << pc_ice << ", " << pc_liq << std::endl;
  std::cout << "  sats = " << sats[0] << ", " << sats[1] << ", " << sats[2] << std::endl;

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
  std::cout << "  dsat_dliq (FD) = " << sats2[0] << ", " << sats2[1] << ", " << sats2[2]
            << std::endl;
  std::cout << "  dsat_dliq = " << sats[0] << ", " << sats[1] << ", " << sats[2] << std::endl;

  p1.saturations(pc_liq, pc_ice + eps, sats2);
  p1.saturations(pc_liq, pc_ice - eps, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2. * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2. * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2. * eps);
  p1.dsaturations_dpc_ice(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3 + 1.e-10);
  std::cout << "  dsat_dice (FD) = " << sats2[0] << ", " << sats2[1] << ", " << sats2[2]
            << std::endl;
  std::cout << "  dsat_dice = " << sats[0] << ", " << sats[1] << ", " << sats[2] << std::endl;

  // unsaturated, barely below freezing
  // -- value
  pc_ice = pcice.CapillaryPressure(273.1, rho);
  pc_liq = p_atm - (-200000);
  p1.saturations(pc_liq, pc_ice, sats);

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

  p1.saturations(pc_liq, pc_ice + eps, sats2);
  p1.saturations(pc_liq, pc_ice - eps, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2 * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2 * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2 * eps);
  p1.dsaturations_dpc_ice(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3 + 1.e-10);


  // unsaturated, well below freezing
  // -- value
  pc_ice = pcice.CapillaryPressure(260., rho);
  pc_liq = p_atm - (-200000);
  p1.saturations(pc_liq, pc_ice, sats);

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

  p1.saturations(pc_liq, pc_ice + eps, sats2);
  p1.saturations(pc_liq, pc_ice - eps, sats3);
  sats2[0] = (sats2[0] - sats3[0]) / (2 * eps);
  sats2[1] = (sats2[1] - sats3[1]) / (2 * eps);
  sats2[2] = (sats2[2] - sats3[2]) / (2 * eps);
  p1.dsaturations_dpc_ice(pc_liq, pc_ice, sats);
  CHECK_CLOSE(sats[0], sats2[0], std::abs(sats[0]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[1], sats2[1], std::abs(sats[1]) / 1.e3 + 1.e-10);
  CHECK_CLOSE(sats[2], sats2[2], std::abs(sats[2]) / 1.e3 + 1.e-10);
}
