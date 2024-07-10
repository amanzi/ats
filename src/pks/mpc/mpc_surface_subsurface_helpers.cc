/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "mpc_surface_subsurface_helpers.hh"
#include "errors.hh"

namespace Amanzi {
namespace MPCHelpers {

void
copySurfaceToSubsurface(const CompositeVector& surf, CompositeVector& sub)
{
  auto parents = surf.getMesh()->getEntityParents<MemSpace_kind::DEVICE>(AmanziMesh::Entity_kind::CELL);
  auto surf_c = surf.viewComponent("cell", false);
  if (sub.hasComponent("face")) {
    auto sub_f = sub.viewComponent("face", false);
    Kokkos::parallel_for(
      "MPCHelpers::copySubsurfaceToSurface", surf_c.extent(0), KOKKOS_LAMBDA(const int& sc) {
        sub_f(parents(sc), 0) = surf_c(sc, 0);
      });

  } else if (sub.hasComponent("boundary_face")) {
    AMANZI_ASSERT(false); // not yet implemented, but could be
  }
}

void
copySubsurfaceToSurface(const CompositeVector& sub, CompositeVector& surf)
{
  auto parents = surf.getMesh()->getEntityParents<MemSpace_kind::DEVICE>(AmanziMesh::Entity_kind::CELL);
  auto surf_c = surf.viewComponent("cell", false);
  if (sub.hasComponent("face")) {
    auto sub_f = sub.viewComponent("face", false);
    Kokkos::parallel_for(
      "MPCHelpers::copySubsurfaceToSurface", surf_c.extent(0), KOKKOS_LAMBDA(const int& sc) {
        surf_c(sc, 0) = sub_f(parents(sc), 0);
      });
  } else if (sub.hasComponent("boundary_face")) {
    AMANZI_ASSERT(false); // not yet implemented, but could be
  }
}

void
mergeSubsurfaceAndSurfacePressure(const CompositeVector& h_prev,
                                  CompositeVector& sub_p,
                                  CompositeVector& surf_p)
{
  double p_atm = 101325.; // get from state?

  auto parents = surf_p.getMesh()->getEntityParents<MemSpace_kind::DEVICE>(AmanziMesh::Entity_kind::CELL);
  auto surf_p_c = surf_p.viewComponent("cell", false);
  auto h_c = h_prev.viewComponent("cell", false);
  if (sub_p.hasComponent("face")) {
    auto sub_p_f = sub_p.viewComponent("face", false);
    Kokkos::parallel_for(
      "MPCHelpers::copySubsurfaceToSurface", surf_p_c.extent(0), KOKKOS_LAMBDA(const int& sc) {
        if (h_c(sc, 0) > 0 && surf_p_c(sc, 0) > p_atm) {
          sub_p_f(parents(sc), 0) = surf_p_c(sc, 0);
        } else {
          surf_p_c(sc, 0) = sub_p_f(parents(sc), 0);
        }
      });
  } else if (sub_p.hasComponent("boundary_face")) {
    AMANZI_ASSERT(false); // not yet implemented, but could be
  }
}

} // namespace MPCHelpers
} // namespace Amanzi
