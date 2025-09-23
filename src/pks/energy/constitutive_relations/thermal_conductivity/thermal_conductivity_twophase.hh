/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Base class of a two-phase Thermal Conductivity relation.

UNITS: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_TC_TWOPHASE_HH_
#define PK_ENERGY_RELATIONS_TC_TWOPHASE_HH_

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class ThermalConductivityTwoPhase {
 public:
  virtual ~ThermalConductivityTwoPhase() {}
  virtual double ThermalConductivity(double porosity, double sat_liq) = 0;
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
