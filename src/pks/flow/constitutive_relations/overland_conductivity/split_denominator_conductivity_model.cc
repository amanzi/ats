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

#include "split_denominator_conductivity_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

SplitDenominatorConductivityModel::SplitDenominatorConductivityModel(Teuchos::ParameterList& plist)
  : plist_(plist)
{
  manning_exp_ = plist_.get<double>("Manning exponent");
}

double
SplitDenominatorConductivityModel::Conductivity(double depth, double slope, double coef)
{
  if (depth <= 0.) return 0.;
  double exponent = manning_exp_ + 1.0;
  return std::pow(std::max(depth, 0.), exponent);
}

double
SplitDenominatorConductivityModel::DConductivityDDepth(double depth, double slope, double coef)
{
  if (depth <= 0.) return 0.;
  double exponent = manning_exp_ + 1.0;
  return std::pow(std::max(depth, 0.), exponent - 1.) * exponent;
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
