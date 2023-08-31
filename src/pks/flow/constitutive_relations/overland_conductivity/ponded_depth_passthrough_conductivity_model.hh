/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth using Manning's model. The denominator in the model is evaluated separately.

*/

#ifndef AMANZI_FLOWRELATIONS_PONDED_DEPTH_PASSTHROUGH_CONDUCTIVITY_MODEL_
#define AMANZI_FLOWRELATIONS_PONDED_DEPTH_PASSTHROUGH_CONDUCTIVITY_MODEL_

#include "Teuchos_ParameterList.hpp"
#include "overland_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

class PondedDepthPassthroughConductivityModel : public OverlandConductivityModel {
 public:
  explicit PondedDepthPassthroughConductivityModel(Teuchos::ParameterList& plist);

  virtual double Conductivity(double depth, double slope, double coef);

  virtual double DConductivityDDepth(double depth, double slope, double coef);

 protected:
  Teuchos::ParameterList plist_;
};

} // namespace Flow
} // namespace Amanzi


#endif
