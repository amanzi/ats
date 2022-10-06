/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! WRMVanGenuchten : water retention model using van Genuchten's parameterization

/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

*/

#ifndef ATS_FLOWRELATIONS_WRM_MACROPORE_
#define ATS_FLOWRELATIONS_WRM_MACROPORE_

#include "Teuchos_ParameterList.hpp"
#include "Spline.hh"

#include "wrm.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Flow {

class WRMMacropore : public WRM {

public:
  explicit WRMMacropore(Teuchos::ParameterList& plist);

  // required methods from the base class
  double k_relative(double saturation);
  double d_k_relative(double saturation);
  double saturation(double pc);
  double d_saturation(double pc);
  double capillaryPressure(double saturation);
  double d_capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }
  //double suction_head(double saturation);
  //double d_suction_head(double saturation);

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList& plist_;

  double a_;  // Macropore parameters: a, dP, Pwe
  double dP_;
  double Pwe_;
  double sr_;  // Macropore residual saturation

  // int function_;
  // double s0_;  // regularization threshold in saturation
  // Amanzi::Utils::Spline fit_kr_;

  // double pc0_;
  // Amanzi::Utils::Spline fit_s_;


  static Utils::RegisteredFactory<WRM,WRMMacropore> factory_;
};

} //namespace
} //namespace

#endif
