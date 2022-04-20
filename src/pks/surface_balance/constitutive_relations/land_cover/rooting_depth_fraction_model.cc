/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Provides a depth-based profile of root density.

#include "rooting_depth_fraction_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
RootingDepthFractionModel::RootingDepthFractionModel(const LandCover& lc)
  : lc_(lc) {}

// main method
double
RootingDepthFractionModel::RootingDepthFraction(double z) const
{
  return (z > lc_.rooting_depth_max) ? 0 :
    0.5 * (lc_.rooting_profile_alpha * exp(-lc_.rooting_profile_alpha*z) + lc_.rooting_profile_beta * exp(-lc_.rooting_profile_beta*z));

}

double
RootingDepthFractionModel::DRootingDepthFractionDDepth(double z) const
{
  return (z > lc_.rooting_depth_max) ? 0 :
    -0.5*pow(lc_.rooting_profile_alpha, 2)*exp(-lc_.rooting_profile_alpha*z) - 0.5*pow(lc_.rooting_profile_beta, 2)*exp(-lc_.rooting_profile_beta*z);
}

} //namespace
} //namespace
} //namespace
