/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

//! WRMPlantChristoffersen : water retention model using van Genuchten's parameterization
/*!
  <ul>Native Spec Example</>
    <ParameterList name="moss" type="ParameterList">
      <Parameter name="region" type="string" value="moss" />
      <Parameter name="wrm type" type="string" value="van Genuchten" />
      <Parameter name="van Genuchten alpha" type="double" value="0.002" />
      <Parameter name="van Genuchten m" type="double" value="0.2" />
      <Parameter name="residual saturation" type="double" value="0.0" />
      <Parameter name="smoothing interval width" type="double" value="100.0" />
    </ParameterList>

*/

#ifndef ATS_FLOWRELATIONS_WRM_PLANTS_CHRISTOFFERSEN_
#define ATS_FLOWRELATIONS_WRM_PLANTS_CHRISTOFFERSEN_

#include "Teuchos_ParameterList.hpp"
#include "wrm.hh"
#include "Factory.hh"
#include <cmath>

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class WRMPlantChristoffersen : public WRM {
 public:
  explicit WRMPlantChristoffersen(Teuchos::ParameterList& plist);

  // required methods from the base class

  double k_relative(double pc);
  double d_k_relative(double pc);
  double saturation(double pc);
  double d_saturation(double pc);
  double capillaryPressure(double saturation);
  double d_capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }

  double operator()(const double& s) const
  {
    return potentialLinear(s) - (potentialSol(s) + potentialP(s));
  }

 private:
  void InitializeFromPlist_();

  /*
  double m_;  // van Genuchten parameters: m, n, alpha
  double n_;
  double l_;
  double alpha_;
  double sr_;  // van Genuchten residual saturation
  */

  double potential(double s) const;
  double potentialLinear(double s) const;
  double potentialSol(double s) const;
  double potentialP(double s) const;

  double d_potential(double s) const;
  double d_potentialLinear(double s) const;
  double d_potentialSol(double s) const;
  double d_potentialP(double s) const;

  Teuchos::ParameterList& plist_;

  double sr_;
  double stlp_;
  double eps_;
  double psi0_;
  double pi0_;
  double psicap_;
  //double star;
  //star = 1.0 - (abs(pi0_)/eps_);
  double sft_;
  double scap_;
  double psitlp_;
  double mcap_;
  double mcapstar_;
  double scapfttrans_;
  double psicapfttrans_;


  struct Tol_ {
    Tol_(double eps)
      : eps_(eps)
    {}
    bool operator()(const double& a, const double& b) const { return std::abs(a - b) <= eps_; }
    double eps_;
  };
  int function_;
  double pc0_;                  // regularization threshold (usually 0 to 500 Pa)
  double a_, b_, factor_dSdPc_; // frequently used constant

  static Utils::RegisteredFactory<WRM, WRMPlantChristoffersen> factory_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
