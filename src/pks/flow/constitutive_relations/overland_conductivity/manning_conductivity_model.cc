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
  manning_exp_ = plist.get<double>("Manning exponent") + 1.0;
}

double ManningConductivityModel::Conductivity(double depth, double slope, double coef) {
  if (depth <= 0.) return 0.;
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  return std::pow(depth, manning_exp_) / scaling;
}

double ManningConductivityModel::DConductivityDDepth(double depth, double slope, double coef) {
  if (depth <= 0.) return 0.;
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  return std::pow(std::max(depth,0.), manning_exp_ - 1.) * manning_exp_ / scaling;
}

} // namespace
} // namespace
