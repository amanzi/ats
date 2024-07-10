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

#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

ManningConductivityModel::ManningConductivityModel(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
{
  slope_regularization_ = plist->get<double>("slope regularization epsilon", 1.e-8);
  manning_exp_ = plist->get<double>("Manning exponent");
  depth_max_ = plist->get<double>("maximum ponded depth [m]", 1.e8);
}


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
