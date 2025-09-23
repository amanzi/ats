/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include "UnitTest++.h"

#include "pc_ice_water.hh"


TEST(pcIL_Derivs)
{
  using namespace Amanzi::Flow::ATS_Physics::Flow;
  double eps = 1.e-4;

  Teuchos::ParameterList plist3;
  PCIceWater pcice(plist3);

  // ensure the above freezing case is 0
  // -- value
  double rho = 1000.;
  double T = 274.;

  CHECK_CLOSE(0., pcice.CapillaryPressure(T, rho), 1.e-10);
  CHECK_CLOSE(0., pcice.DCapillaryPressureDT(T, rho), 1.e-10);
  CHECK_CLOSE(0., pcice.DCapillaryPressureDRho(T, rho), 1.e-10);

  T = 270.;
  double dpdT =
    (pcice.CapillaryPressure(T + eps, rho) - pcice.CapillaryPressure(T - eps, rho)) / (2 * eps);
  double dpdrho =
    (pcice.CapillaryPressure(T, rho + eps) - pcice.CapillaryPressure(T, rho - eps)) / (2 * eps);

  double pc = pcice.CapillaryPressure(T, rho);
  CHECK_CLOSE(dpdT, pcice.DCapillaryPressureDT(T, rho), pc * 1.e-6);
  CHECK_CLOSE(dpdrho, pcice.DCapillaryPressureDRho(T, rho), pc * 1.e-6);
}
