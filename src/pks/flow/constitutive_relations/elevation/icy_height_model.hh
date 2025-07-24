/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_ICY_HEIGHT_MODEL_
#define AMANZI_FLOWRELATIONS_ICY_HEIGHT_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class IcyHeightModel {
 public:
  explicit IcyHeightModel(Teuchos::ParameterList& plist)
    : plist_(plist)
  {}

  double Height(double pres, double eta, double rhol, double rhoi, double p_atm, double g_z)
  {
    return (pres - p_atm) / ((eta * rhol + (1. - eta) * rhoi) * g_z);
  }

  double DHeightDPressure(double pres,
                          double eta,
                          double rhol,
                          double rhoi,
                          double p_atm,
                          double g_z)
  {
    return 1. / ((eta * rhol + (1. - eta) * rhoi) * g_z);
  }

  double DHeightDRho_l(double pres, double eta, double rhol, double rhoi, double p_atm, double g_z)
  {
    double denom = ((eta * rhol + (1.0 - eta) * rhoi) * g_z);
    return -(pres - p_atm) / (denom * denom) * g_z * eta;
  }

  double DHeightDRho_i(double pres, double eta, double rhol, double rhoi, double p_atm, double g_z)
  {
    double denom = ((eta * rhol + (1.0 - eta) * rhoi) * g_z);
    return -(pres - p_atm) / (denom * denom) * g_z * (1. - eta);
  }

  double DHeightDEta(double pres, double eta, double rhol, double rhoi, double p_atm, double g_z)
  {
    double denom = ((eta * rhol + (1.0 - eta) * rhoi) * g_z);
    return -(pres - p_atm) / (denom * denom) * g_z * (rhol - rhoi);
  }

 protected:
  Teuchos::ParameterList plist_;
};

} // namespace Flow
} // namespace Amanzi

#endif
