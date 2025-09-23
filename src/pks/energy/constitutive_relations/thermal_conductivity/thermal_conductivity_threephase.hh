/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Base class of a three-phase Thermal Conductivity relation.

UNITS: ????
------------------------------------------------------------------------- */

#ifndef PK_ENERGY_RELATIONS_TC_THREEPHASE_HH_
#define PK_ENERGY_RELATIONS_TC_THREEPHASE_HH_

#include "dbc.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class ThermalConductivityThreePhase {
 public:
  virtual ~ThermalConductivityThreePhase() {}

  virtual double ThermalConductivity(double porosity,
                                     double sat_liq,
                                     double sat_ice,
                                     double temp) = 0;
  virtual double DThermalConductivity_DPorosity(double porosity,
                                                double sat_liq,
                                                double sat_ice,
                                                double temp)
  {
    AMANZI_ASSERT(false);
    return 0.;
  }
  virtual double DThermalConductivity_DSaturationLiquid(double porosity,
                                                        double sat_liq,
                                                        double sat_ice,
                                                        double temp)
  {
    AMANZI_ASSERT(false);
    return 0.;
  }
  virtual double DThermalConductivity_DSaturationIce(double porosity,
                                                     double sat_liq,
                                                     double sat_ice,
                                                     double temp)
  {
    AMANZI_ASSERT(false);
    return 0.;
  }
  virtual double DThermalConductivity_DTemperature(double porosity,
                                                   double sat_liq,
                                                   double sat_ice,
                                                   double temp)
  {
    AMANZI_ASSERT(false);
    return 0.;
  }
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
