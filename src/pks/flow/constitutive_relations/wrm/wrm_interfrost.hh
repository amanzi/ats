/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Amanzi needs to fix its WRM interface... the current lack of factory or
  common interface is stupid.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef FLOWRELATIONS_WRM_INTERFROST_
#define FLOWRELATIONS_WRM_INTERFROST_

#include "Teuchos_ParameterList.hpp"
#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
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
} // namespace Amanzi

#endif
