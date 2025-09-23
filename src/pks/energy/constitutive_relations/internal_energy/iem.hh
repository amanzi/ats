/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Internal energy model -- function of temperature only.

UNITS: J/{mol/kg}
------------------------------------------------------------------------- */

#ifndef AMANZI_ENERGYRELATIONS_IEM_
#define AMANZI_ENERGYRELATIONS_IEM_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

class IEM {
 public:
  virtual ~IEM() {}

  // IEM(Teuchos::ParameterList& plist);
  virtual bool IsMolarBasis() = 0;
  virtual double InternalEnergy(double temp) = 0;
  virtual double DInternalEnergyDT(double temp) = 0;
};

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi

#endif
