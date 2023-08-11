/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#ifndef PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_

#include "CompositeVector.hh"

namespace Amanzi {
namespace MPCHelpers {

void
copySurfaceToSubsurface(const CompositeVector& surf, CompositeVector& sub);

void
copySubsurfaceToSurface(const CompositeVector& sub, CompositeVector& surf);

void
mergeSubsurfaceAndSurfacePressure(const CompositeVector& kr_surf,
                                  CompositeVector& sub_p,
                                  CompositeVector& surf_p);

} // namespace MPCHelpers
} // namespace Amanzi


#endif
