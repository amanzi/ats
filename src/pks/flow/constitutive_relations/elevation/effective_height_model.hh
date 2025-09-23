/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates height(pressure)

*/

#ifndef AMANZI_FLOWRELATIONS_EFFECTIVE_HEIGHT_MODEL_
#define AMANZI_FLOWRELATIONS_EFFECTIVE_HEIGHT_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class EffectiveHeightModel {
 public:
  explicit EffectiveHeightModel(Teuchos::ParameterList& plist)
    : plist_(plist)
  {
    InitializeFromPList_();
  }

  double EffectiveHeight(double h)
  {
    if (h <= 0.) {
      return 0.;
    } else if (h >= smoothing_width_) {
      return h - 0.5 * smoothing_width_;
    } else {
      return 0.5 * h - 0.5 * std::sin(h / smoothing_width_ * pi_) * smoothing_width_ / pi_;
    }
  }

  double DEffectiveHeightDHeight(double h)
  {
    if (h <= 0.) {
      return 0.;
    } else if (h >= smoothing_width_) {
      return 1.;
    } else {
      return 0.5 - 0.5 * std::cos(h / smoothing_width_ * pi_);
    }
  }

 protected:
  void InitializeFromPList_()
  {
    smoothing_width_ = plist_.get<double>("smoothing width [m]", 0.01);
    pi_ = 4 * std::atan(1.0);
  }

  Teuchos::ParameterList plist_;
  double smoothing_width_;
  double pi_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
