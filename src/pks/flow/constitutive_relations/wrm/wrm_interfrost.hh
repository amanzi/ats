/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

This model was used as a part of the INTERFROST code comparison effort which
led to the paper by Grenier et al 2018 AWR.  It is not intended to be used in
most ATS runs.

*/
#ifndef FLOWRELATIONS_WRM_INTERFROST_
#define FLOWRELATIONS_WRM_INTERFROST_

#include "Teuchos_ParameterList.hpp"
#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class WRMInterfrost : public WRM {
 public:
  explicit WRMInterfrost(Teuchos::ParameterList& plist) {}

  // required methods from the base class
  double k_relative(double sat) { return std::pow(10, -50 * 0.37 * (1 - sat)); }
  double d_k_relative(double pc) { return 0.; }
  double saturation(double pc)
  {
    AMANZI_ASSERT(0);
    return 0.;
  }
  double d_saturation(double pc)
  {
    AMANZI_ASSERT(0);
    return 0.;
  }
  double capillaryPressure(double saturation) { return saturation; }
  double d_capillaryPressure(double saturation) { return 1.; }
  double residualSaturation()
  {
    AMANZI_ASSERT(0);
    return 0.;
  }

 private:
  static Utils::RegisteredFactory<WRM, WRMInterfrost> factory_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
