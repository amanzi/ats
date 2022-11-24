/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "test_snow_dist.hh"

namespace Amanzi {

RegisteredPKFactory_ATS<TestSnowDist> TestSnowDist::reg_("snow distribution test");

} // namespace Amanzi
