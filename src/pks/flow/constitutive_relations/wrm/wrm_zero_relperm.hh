/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  WRM which calls another WRM for saturation but sets 0 rel perm.

*/

#ifndef AMANZI_FLOWRELATIONS_WRM_ZERO_RELPERM_HH_
#define AMANZI_FLOWRELATIONS_WRM_ZERO_RELPERM_HH_

#include "Teuchos_ParameterList.hpp"
#include "wrm.hh"
#include "Factory.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {

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

  static Utils::RegisteredFactory<WRM, WRMZeroRelPerm> factory_;
};

} // namespace Flow
} // namespace Amanzi

#endif
