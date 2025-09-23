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

#include "ponded_depth_passthrough_conductivity_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

PondedDepthPassthroughConductivityModel::PondedDepthPassthroughConductivityModel(
  Teuchos::ParameterList& plist)
  : plist_(plist)
{}

double
PondedDepthPassthroughConductivityModel::Conductivity(double depth, double slope, double coef)
{
  return depth;
}

double
PondedDepthPassthroughConductivityModel::DConductivityDDepth(double depth,
                                                             double slope,
                                                             double coef)
{
  return 1;
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
