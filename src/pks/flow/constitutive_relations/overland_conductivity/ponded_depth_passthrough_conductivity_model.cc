/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
 Evaluates the conductivity of surface flow as a function of ponded
 depth using Manning's model. The denominator in the model is evaluated separately.
 
*/

#include "ponded_depth_passthrough_conductivity_model.hh"

namespace Amanzi {
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
} // namespace Amanzi
