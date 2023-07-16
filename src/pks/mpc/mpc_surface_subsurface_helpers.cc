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

void
CopySurfaceToSubsurface(const CompositeVector& surf, CompositeVector& sub)
{
  const Epetra_MultiVector& surf_c = *surf.ViewComponent("cell", false);

  if (sub.HasComponent("face")) {
    DomainFaceSetter sub_f(*sub.Mesh(), *sub.ViewComponent("face", false));

    for (int sc = 0; sc != surf_c.MyLength(); ++sc) {
      AmanziMesh::Entity_ID f = surf.Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
      sub_f.set<AmanziMesh::FACE>(f, surf_c[0][sc]);
    }

  } else if (sub.HasComponent("boundary_face")) {
    DomainFaceSetter sub_f(*sub.Mesh(), *sub.ViewComponent("boundary_face", false));

    for (int sc = 0; sc != surf_c.MyLength(); ++sc) {
      AmanziMesh::Entity_ID f = surf.Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
      sub_f.set<AmanziMesh::BOUNDARY_FACE>(f, surf_c[0][sc]);
    }
  }
}

void
CopySubsurfaceToSurface(const CompositeVector& sub, CompositeVector& surf)
{
  //  const Epetra_MultiVector& sub_f = *sub.ViewComponent("face",false);
  Epetra_MultiVector& surf_c = *surf.ViewComponent("cell", false);

  if (sub.HasComponent("face")) {
    DomainFaceGetter sub_f(*sub.Mesh(), *sub.ViewComponent("face", false));

    for (int sc = 0; sc != surf_c.MyLength(); ++sc) {
      AmanziMesh::Entity_ID f = surf.Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
      surf_c[0][sc] = sub_f.get<AmanziMesh::FACE>(f);
    }

  } else if (sub.HasComponent("boundary_face")) {
    DomainFaceGetter sub_f(*sub.Mesh(), *sub.ViewComponent("boundary_face", false));

    for (int sc = 0; sc != surf_c.MyLength(); ++sc) {
      AmanziMesh::Entity_ID f = surf.Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
      surf_c[0][sc] = sub_f.get<AmanziMesh::BOUNDARY_FACE>(f);
    }
  }
}

void
MergeSubsurfaceAndSurfacePressure(const CompositeVector& h_prev,
                                  CompositeVector& sub_p,
                                  CompositeVector& surf_p)
{
  Epetra_MultiVector& surf_p_c = *surf_p.ViewComponent("cell", false);
  const Epetra_MultiVector& h_c = *h_prev.ViewComponent("cell", false);
  double p_atm = 101325.;

  if (sub_p.HasComponent("face")) {
    DomainFaceSetter sub_p_f(*sub_p.Mesh(), *sub_p.ViewComponent("face", false));

    for (unsigned int sc = 0; sc != surf_p_c.MyLength(); ++sc) {
      AmanziMesh::Entity_ID f = surf_p.Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
      if (h_c[0][sc] > 0. && surf_p_c[0][sc] > p_atm) {
        sub_p_f.set<AmanziMesh::FACE>(f, surf_p_c[0][sc]);
      } else {
        surf_p_c[0][sc] = sub_p_f.get<AmanziMesh::FACE>(f);
      }
    }

  } else if (sub_p.HasComponent("boundary_face")) {
    DomainFaceSetter sub_p_f(*sub_p.Mesh(), *sub_p.ViewComponent("boundary_face", false));

    for (unsigned int sc = 0; sc != surf_p_c.MyLength(); ++sc) {
      AmanziMesh::Entity_ID f = surf_p.Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
      if (h_c[0][sc] > 0. && surf_p_c[0][sc] > p_atm) {
        sub_p_f.set<AmanziMesh::BOUNDARY_FACE>(f, surf_p_c[0][sc]);
      } else {
        surf_p_c[0][sc] = sub_p_f.get<AmanziMesh::BOUNDARY_FACE>(f);
      }
    }
  }
}


} // namespace Amanzi
