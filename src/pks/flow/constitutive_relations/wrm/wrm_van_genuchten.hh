/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! WRMVanGenuchten : water retention model using van Genuchten's parameterization

/*!
  <ul>Native Spec Example</>
    <ParameterList name="moss" type="ParameterList">
      <Parameter name="region" type="string" value="moss" />
      <Parameter name="WRM Type" type="string" value="van Genuchten" />
      <Parameter name="van Genuchten alpha" type="double" value="0.002" />
      <Parameter name="van Genuchten m" type="double" value="0.2" />
      <Parameter name="residual saturation" type="double" value="0.0" />
      <Parameter name="smoothing interval width [saturation]" type="double" value=".05" />
    </ParameterList>

*/

#ifndef ATS_FLOWRELATIONS_WRM_VAN_GENUCHTEN_
#define ATS_FLOWRELATIONS_WRM_VAN_GENUCHTEN_

#include "Teuchos_ParameterList.hpp"
#include "Spline.hh"

#include "wrm.hh"
#include "factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMVanGenuchten : public WRM {

public:
  explicit WRMVanGenuchten(Teuchos::ParameterList& plist);

  // required methods from the base class
  double k_relative(double saturation);
  double d_k_relative(double saturation);
  double saturation(double pc);
  double d_saturation(double pc);
  double capillaryPressure(double saturation);
  double d_capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList& plist_;

  double m_;  // van Genuchten parameters: m, n, alpha
  double n_;
  double l_;
  double alpha_;
  double sr_;  // van Genuchten residual saturation

  int function_;
  double s0_;  // regularization threshold in saturation
  Amanzi::Utils::Spline fit_;
  

  static Utils::RegisteredFactory<WRM,WRMVanGenuchten> factory_;
};

} //namespace
} //namespace
} //namespace

#endif
