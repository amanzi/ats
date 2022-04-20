/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Provides a depth-based profile of root density.
/*!

Sets the root fraction as a function of depth,

.. math:
   F_root =  ( \alpha \; exp(-\alpha z) + \beta \; exp(-\beta z) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.

Note that all three parameters, a, b, and the cutoff, are provided in the
LandCover type.

*/

#pragma once

#include "LandCover.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class RootingDepthFractionModel {

 public:
  explicit
  RootingDepthFractionModel(const LandCover& lc);

  double RootingDepthFraction(double z) const;

  double DRootingDepthFractionDDepth(double z) const;

 protected:
  const LandCover& lc_;

};

} //namespace
} //namespace
} //namespace

