/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  A capillary pressure model based upon something other than p_atm - p.

*/

#ifndef AMANZI_FLOW_RELATIONS_PC_LIQ_ATM_
#define AMANZI_FLOW_RELATIONS_PC_LIQ_ATM_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class PCLiqAtm {
 public:
  explicit PCLiqAtm(Teuchos::ParameterList& plist) {}

  // required methods from the base class
  double CapillaryPressure(double p, double p_atm) { return p_atm - p; }
  double DCapillaryPressureDp(double p, double p_atm) { return -1.; }
  double DCapillaryPressureDpatm(double p, double p_atm) { return 1.; }
};

} // namespace Flow
} // namespace Amanzi

#endif
