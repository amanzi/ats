/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Helper functions for the subgrid topography model.
//! Note these evaluate equation 7 of Jan et al WRR 2018, and its derivative.

#pragma once

namespace Amanzi {
namespace Flow {

double subgrid_VolumetricDepth(double depth, double del_max, double del_ex);
double subgrid_DVolumetricDepth_DDepth(double depth, double del_max, double del_ex);

} // namespace
} // namespace
