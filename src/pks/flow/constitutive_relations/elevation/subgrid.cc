/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

#include <cmath>
#include "subgrid.hh"

namespace Amanzi {
namespace Flow {

double
subgrid_VolumetricDepth(double depth, double del_max, double del_ex)
{
  return depth >= del_max ? depth - del_ex:
      std::pow(depth/del_max, 2) * (2*del_max - 3*del_ex)
      + std::pow(depth/del_max,3) * (2*del_ex - del_max);
}

double
subgrid_DVolumetricDepth_DDepth(double depth, double del_max, double del_ex)
{
  return depth >= del_max ? 1 :
      2 * depth/del_max * (2*del_max - 3*del_ex) / del_max
      + 3 * std::pow(depth/del_max,2) * (2*del_ex - del_max) / del_max;
}  

} // namespace
} // namespace
