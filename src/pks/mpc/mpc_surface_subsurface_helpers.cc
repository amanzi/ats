#include "mpc_surface_subsurface_helpers.hh"
#include "errors.hh"

namespace Amanzi {

void
CopySurfaceToSubsurface(const CompositeVector& surf,
                        const Teuchos::Ptr<CompositeVector>& sub)
{
  const Epetra_MultiVector& surf_c = *surf.ViewComponent("cell",false);

  for (unsigned int sc=0; sc!=surf_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
        surf.Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
    SetDomainFaceValue(*sub, f, surf_c[0][sc]);
  }
}

void
CopySubsurfaceToSurface(const CompositeVector& sub,
                        const Teuchos::Ptr<CompositeVector>& surf)
{
  //  const Epetra_MultiVector& sub_f = *sub.ViewComponent("face",false);
  Epetra_MultiVector& surf_c = *surf->ViewComponent("cell",false);

  for (unsigned int sc=0; sc!=surf_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
        surf->Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
    surf_c[0][sc] = GetDomainFaceValue(sub, f);
  }
}

void
MergeSubsurfaceAndSurfacePressure(const CompositeVector& h_prev,
				  const Teuchos::Ptr<CompositeVector>& sub_p,
				  const Teuchos::Ptr<CompositeVector>& surf_p)
{
  Epetra_MultiVector& surf_p_c = *surf_p->ViewComponent("cell",false);
  const Epetra_MultiVector& h_c = *h_prev.ViewComponent("cell",false);
  double p_atm = 101325.;

  for (unsigned int sc=0; sc!=surf_p_c.MyLength(); ++sc) {
    AmanziMesh::Entity_ID f =
        surf_p->Mesh()->entity_get_parent(AmanziMesh::CELL, sc);
    if (h_c[0][sc] > 0. && surf_p_c[0][sc] > p_atm) {
      SetDomainFaceValue(*sub_p, f,  surf_p_c[0][sc]);
    } else {
      surf_p_c[0][sc] = GetDomainFaceValue(*sub_p, f);
    }
  }
}

double
GetDomainFaceValue(const CompositeVector& sub_p, int f)
{
  std::string face_entity;
  if (sub_p.HasComponent("face")) {
    face_entity = "face";
  } else if (sub_p.HasComponent("boundary_face")) {
    face_entity = "boundary_face";
  } else {
    Errors::Message message("Subsurface vector does not have face component.");
    Exceptions::amanzi_throw(message);
  }

  if (face_entity == "face") {
    const Epetra_MultiVector& vec = *sub_p.ViewComponent(face_entity, false);
    return vec[0][f];;
  } else if (face_entity == "boundary_face") {
    int bf = sub_p.Mesh()->exterior_face_map(false).LID(sub_p.Mesh()->face_map(false).GID(f));
    const Epetra_MultiVector& vec = *sub_p.ViewComponent(face_entity, false);
    return vec[0][bf];
  } else {
    return std::numeric_limits<double>::quiet_NaN();
  }
}

void
SetDomainFaceValue(CompositeVector& sub_p, int f, double value)
{
  std::string face_entity;
  if (sub_p.HasComponent("face")) {
    face_entity = "face";
  } else if (sub_p.HasComponent("boundary_face")) {
    face_entity = "boundary_face";
  } else {
    Errors::Message message("Subsurface vector does not have face component.");
    Exceptions::amanzi_throw(message);
  }

  if (face_entity == "face") {
    Epetra_MultiVector& vec = *sub_p.ViewComponent(face_entity, false);
    vec[0][f] = value;
  } else if (face_entity == "boundary_face") {
    int bf = sub_p.Mesh()->exterior_face_map(false).LID(sub_p.Mesh()->face_map(false).GID(f));
    Epetra_MultiVector& vec = *sub_p.ViewComponent(face_entity, false);
    vec[0][bf] = value;
  }

}


} // namespace
