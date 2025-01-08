/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth and surface slope using Manning's model.

*/

#ifndef AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_
#define AMANZI_FLOWRELATIONS_MANNING_CONDUCTIVITY_MODEL_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Flow {

class ManningConductivityModel {
 public:
  explicit ManningConductivityModel(Teuchos::ParameterList& plist);

  double Conductivity(double depth, double slope, double coef);
  double DConductivityDDepth(double depth, double slope, double coef);

 protected:
  double slope_regularization_;
  double manning_exp_;
  double depth_max_;
};

} // namespace Flow
} // namespace Amanzi

#endif
