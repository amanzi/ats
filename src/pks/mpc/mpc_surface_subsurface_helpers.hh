#ifndef PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_
#define PKS_MPC_SURFACE_SUBSURFACE_HELPERS_HH_

#include "CompositeVector.hh"

namespace Amanzi {

void
CopySurfaceToSubsurface(const CompositeVector& surf, CompositeVector& sub);

void
CopySubsurfaceToSurface(const CompositeVector& sub, CompositeVector& surf);

void
MergeSubsurfaceAndSurfacePressure(const CompositeVector& kr_surf,
                                  CompositeVector& sub_p,
                                  CompositeVector& surf_p);
double
GetDomainFaceValue(const CompositeVector& sub_p, int f);

void
SetDomainFaceValue(CompositeVector& sub_p, int f, double value);


} // namespace Amanzi


#endif
