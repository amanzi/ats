/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth and surface slope using Manning's model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

ManningConductivityModel::ManningConductivityModel(Teuchos::ParameterList& plist)
{
  slope_regularization_ = plist.get<double>("slope regularization epsilon", 1.e-8);
  manning_exp_ = plist.get<double>("Manning exponent");
  depth_max_ = plist.get<double>("maximum ponded depth [m]", 1.e8);
}

double ManningConductivityModel::Conductivity(double depth, double slope, double coef) {
  if (depth <= 0.) return 0.;
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  return depth * std::pow(std::min(depth, depth_max_), manning_exp_) / scaling;
}

double ManningConductivityModel::DConductivityDDepth(double depth, double slope, double coef) {
  if (depth <= 0.) return 0.;
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  if (depth > depth_max_) {
    return std::pow(depth_max_, manning_exp_) / scaling;
  } else {
    return std::pow(depth, manning_exp_) * (manning_exp_+1) / scaling;
  }
}

} // namespace
} // namespace
