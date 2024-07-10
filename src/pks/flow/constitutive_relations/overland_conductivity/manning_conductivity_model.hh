/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth and surface slope using Manning's model.

*/

#ifndef AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_
#define AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_

#include "Kokkos_Core.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ManningConductivityModel {
 public:
  explicit ManningConductivityModel(const Teuchos::RCP<Teuchos::ParameterList>& plist);

  KOKKOS_INLINE_FUNCTION
  double Conductivity(double depth, double slope, double coef) const {
    if (depth <= 0.) return 0.;
    double scaling = coef * Kokkos::sqrt(Kokkos::max(slope, slope_regularization_));
    return Kokkos::pow(Kokkos::min(depth, depth_max_), manning_exp_) / scaling;
  }

  KOKKOS_INLINE_FUNCTION
  double DConductivityDDepth(double depth, double slope, double coef) const {
    if (depth <= 0.) return 0.;
    double scaling = coef * Kokkos::sqrt(Kokkos::max(slope, slope_regularization_));
    if (depth > depth_max_) {
      return 0.;
    } else {
      return manning_exp_ * Kokkos::pow(depth, manning_exp_ - 1) / scaling;
    }
  }

 protected:
  double slope_regularization_;
  double manning_exp_;
  double depth_max_;
};

} // namespace Relations
} // namespace Flow
} // namespace Amanzi

#endif
