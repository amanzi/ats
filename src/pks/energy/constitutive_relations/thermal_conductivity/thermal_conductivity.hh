/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Base of a Thermal Conductivity relation.

UNITS: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_HH_
#define PK_ENERGY_RELATIONS_THERMAL_CONDUCTIVITY_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class TwophaseThermalConductivity {
 public:
  TwophaseThermalConductivity(Teuchos::ParameterList& plist);

  double ThermalConductivity(double porosity, double sat_liq);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;

  double eps_;
  double alpha_;
  double k_soil_;
  double k_liquid_;
  double k_gas_;
  double d_;
  double rho_soil_;
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
