/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#pragma once

#include <limits>
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

#define AMANZI_NAN_SIGNALING
#ifdef AMANZI_NAN_SIGNALING
static const double NaN = std::numeric_limits<double>::signaling_NaN();
#else
static const double NaN = std::numeric_limits<double>::quiet_NaN();
#endif

} // namespace
} // namespace
} // namespace

