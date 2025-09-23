/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <iostream>
#include "UnitTest++.h"

#include "wrm_plants_christoffersen.hh"

TEST(plantChristoffersen)
{
  using namespace Amanzi::ATS_Physics::Flow;

  double sr = 0.46567;
  double stlp = 0.8975;
  double eps = 15.90459885634534;
  double psi0 = -0.08;
  double pi0 = -1.93116;
  double psicap = -0.39;
  double p_atm = 101325.0;
  double pc;


  Teuchos::ParameterList plist;
  plist.set("residual saturation [-]", sr);
  plist.set("saturation at turgor loss [-]", stlp);
  plist.set("bulk elastic modulus [MPa]", eps);
  plist.set("water potential at full saturation [Pa]", psi0);
  plist.set("osmotic potential at full turgor [Pa]", pi0);
  plist.set("water potential at full turgor [Pa]", psicap);
  plist.set("smoothing interval width [saturation]", 0.0);
  WRMPlantChristoffersen pC(plist);
  /*
  // check k_relative for p = 2*p_atm
  double pc = -p_atm;
  CHECK_EQUAL(pC.k_relative(pc), 1.0);

  // check k_relative for p = 0
  pc = p_atm;
  double se = std::pow(1.0 + std::pow(alpha * pc, 1.0 / (1.0-m)),-m);
  CHECK_CLOSE(vG.k_relative(pc),
              sqrt(se) * std::pow(1.0 - std::pow(1.0 - std::pow(se, 1.0/m), m), 2.0), 1e-15);
  */
  // check saturation for p = 2*p_atm
  pc = p_atm;
  CHECK_CLOSE(pC.saturation(pc), 0.995698885081, 1e-11);

  // check saturation for p = 0
  pc = 100 * p_atm;
  CHECK_CLOSE(pC.saturation(pc), 0.559347194579, 1e-11);

  // check derivative of saturation(p) at p=2*p_atm
  pc = 100 * p_atm;
  //  CHECK_ClOSE(pC.d_saturation(pc),-9.24522028989e-09,1e-14);

  // check derivative of saturation(p) at p=0
  pc = p_atm;
  CHECK_CLOSE(pC.d_saturation(pc), -2.01693548418e-07, 1e-12);

  pc = 100 * p_atm;
  CHECK_CLOSE(pC.d_saturation(pc), -9.24522028989e-09, 1e-15);

  // check capillary pressure at p = 2*p_atm
  pc = p_atm;
  CHECK_CLOSE(pC.capillaryPressure(0.9), 2104501.95778, 1e-4);

  // check capillary pressure at p = 0.
  pc = 100 * p_atm;
  CHECK_CLOSE(pC.capillaryPressure(pC.saturation(pc)), pc, 1e-2);

  // check d_capillaryPressure at p = 0
  pc = p_atm;
  CHECK_CLOSE(pC.d_capillaryPressure(pC.saturation(pc)), -4958016.79327, 1.);
}
