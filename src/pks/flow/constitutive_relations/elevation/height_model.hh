/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates height(pressure)

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_HEIGHT_MODEL_
#define AMANZI_FLOWRELATIONS_HEIGHT_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class HeightModel {
 public:
  explicit HeightModel(Teuchos::ParameterList& plist) : plist_(plist) {}

  double Height(double pres, double rho, double p_atm, double g_z)
  {
    return (pres - p_atm) / (rho * g_z);
  }

  double DHeightDPressure(double pres, double rho, double p_atm, double g_z)
  {
    return 1. / (rho * g_z);
  }

  double DHeightDRho(double pres, double rho, double p_atm, double g_z)
  {
    return -(pres - p_atm) / (rho * rho * g_z);
  }

 protected:
  Teuchos::ParameterList plist_;
};

} // namespace Flow
} // namespace Amanzi

#endif
