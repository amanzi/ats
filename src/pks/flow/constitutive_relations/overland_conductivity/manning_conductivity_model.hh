/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth and surface slope using Manning's model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_
#define AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel {
public:
  explicit
  ManningConductivityModel(Teuchos::ParameterList& plist);

  double Conductivity(double depth, double slope, double coef);
  double DConductivityDDepth(double depth, double slope, double coef);

protected:
  double slope_regularization_;
  double manning_exp_;
  double depth_max_;
};

} // namespace
} // namespace

#endif
