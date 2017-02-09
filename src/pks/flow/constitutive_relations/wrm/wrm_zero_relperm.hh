/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  WRM which calls another WRM for saturation but sets 0 rel perm.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_WRM_ZERO_RELPERM_HH_
#define AMANZI_FLOWRELATIONS_WRM_ZERO_RELPERM_HH_

#include "Teuchos_ParameterList.hpp"
#include "wrm.hh"
#include "factory.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

class WRMZeroRelPerm : public WRM {

 public:
  explicit WRMZeroRelPerm(Teuchos::ParameterList& plist);

  double k_relative(double s) { return 0.0; }
  double d_k_relative(double s) { return 0.0; }
  double saturation(double pc) { return wrm_->saturation(pc); }
  double d_saturation(double pc) { return wrm_->d_saturation(pc); }
  double capillaryPressure(double sat) { return wrm_->capillaryPressure(sat); }
  double d_capillaryPressure(double sat) { return wrm_->d_capillaryPressure(sat); }
  double residualSaturation() { return wrm_->residualSaturation(); }

 private:
  void InitializeFromPlist_();

  Teuchos::ParameterList plist_;
  Teuchos::RCP<WRM> wrm_;

  static Utils::RegisteredFactory<WRM,WRMZeroRelPerm> factory_;
};

} // namespace
} // namespace
} // namespace

#endif

