/* -*-  mode: c++; indent-tabs-mode: nil -*- */

// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
// Author: Ethan Coon (ecoon@lanl.gov)
//
// Scheme for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#include "Tensor.hh"
#include "CompositeVector.hh"
#include "State.hh"
#include "upwind_gravity_flux.hh"

namespace Amanzi {
namespace Operators {

UpwindGravityFlux::UpwindGravityFlux(std::string pkname,
        std::string cell_coef,
        std::string face_coef,
        const Teuchos::RCP<std::vector<WhetStone::Tensor> > K) :
    pkname_(pkname),
    cell_coef_(cell_coef),
    face_coef_(face_coef),
    K_(K) {};

void UpwindGravityFlux::Update(const Teuchos::Ptr<State>& S,
                               const Teuchos::Ptr<Debugger>& db) {

  Teuchos::RCP<const CompositeVector> cell = S->GetFieldData(cell_coef_);
  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");
  Teuchos::RCP<CompositeVector> face = S->GetFieldData(face_coef_, pkname_);
  CalculateCoefficientsOnFaces(*cell, "cell", *g_vec, face.ptr(), "face");
};
  

void UpwindGravityFlux::Update(const Teuchos::Ptr<State>& S,
                               const std::string cell_key,
                               const std::string cell_component,
                               const std::string face_key,
                               const std::string face_component,
                               const Teuchos::Ptr<Debugger>& db) {

  Teuchos::RCP<const CompositeVector> cell = S->GetFieldData(cell_key);
  Teuchos::RCP<const Epetra_Vector> g_vec = S->GetConstantVectorData("gravity");
  Teuchos::RCP<CompositeVector> face = S->GetFieldData(face_key, pkname_);
  CalculateCoefficientsOnFaces(*cell,cell_component,*g_vec, face.ptr(), face_component);

};


void UpwindGravityFlux::CalculateCoefficientsOnFaces(
        const CompositeVector& cell_coef,
        const std::string cell_component,
        const Epetra_Vector& g_vec,
        const Teuchos::Ptr<CompositeVector>& face_coef,
        const std::string face_component) {

  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  double flow_eps = 1.e-10;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh = face_coef->Mesh();

  // set up gravity
  AmanziGeometry::Point gravity(g_vec.MyLength());
  for (int i=0; i!=g_vec.MyLength(); ++i) gravity[i] = g_vec[i];

  // initialize the face coefficients
  face_coef->ViewComponent(face_component,true)->PutScalar(0.0);
  if (face_coef->HasComponent(cell_component)) {
    face_coef->ViewComponent(cell_component,true)->PutScalar(1.0);
  }

  // Note that by scattering, and then looping over all Parallel_type::ALL cells, we
  // end up getting the correct upwind values in all faces (owned or
  // not) bordering an owned cell.  This is the necessary data for
  // making the local matrices in MFD, so there is no need to
  // communicate the resulting face coeficients.

  // communicate ghosted cells
  cell_coef.ScatterMasterToGhosted(cell_component);

  Epetra_MultiVector& face_coef_v = *face_coef->ViewComponent(face_component,true);
  const Epetra_MultiVector& cell_coef_v = *cell_coef.ViewComponent(cell_component,true);


  for (unsigned int c=0; c!=cell_coef.size(cell_component, true); ++c) {
    mesh->cell_get_faces_and_dirs(c, &faces, &dirs);
    AmanziGeometry::Point Kgravity = (*K_)[c] * gravity;

    for (unsigned int n=0; n!=faces.size(); ++n) {
      int f = faces[n];

      const AmanziGeometry::Point& normal = mesh->face_normal(f);
      if ((normal * Kgravity) * dirs[n] >= flow_eps) {
        face_coef_v[0][f] = cell_coef_v[0][c];
      } else if (std::abs((normal * Kgravity) * dirs[n]) < flow_eps) {
        face_coef_v[0][f] += cell_coef_v[0][c] / 2.;
      }
    }
  }
};


} //namespace
} //namespace
