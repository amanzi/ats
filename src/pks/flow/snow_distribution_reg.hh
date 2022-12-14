/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Author: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#include "snow_distribution.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<SnowDistribution> SnowDistribution::reg_("snow distribution");

} // namespace Flow
} // namespace Amanzi
