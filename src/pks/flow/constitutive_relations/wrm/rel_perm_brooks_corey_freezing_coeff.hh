/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Bo Gao (gaob@ornl.gov)
*/

#pragma once
#include <cmath>

namespace Amanzi {
namespace ATS_Physics {          
namespace Flow {
namespace BrooksCoreyFrzCoef {


inline double
frzcoef(double sl, double sg, double omega)
{
  return 1. - std::exp(-omega * (sl + sg)) + std::exp(-omega);
}


inline double
d_frzcoef_dsl(double sl, double sg, double omega)
{
  return omega * std::exp(-omega * (sl + sg));
}


inline double
d_frzcoef_dsg(double sl, double sg, double omega)
{
  return omega * std::exp(-omega * (sl + sg));
}


} // namespace BrooksCoreyFrzCoef
} // namespace Flow
}
} // namespace Amanzi
