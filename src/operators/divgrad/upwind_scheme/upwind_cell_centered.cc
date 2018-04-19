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

#include "Mesh.hh"

#include "CompositeVector.hh"
#include "State.hh"
#include "upwind_cell_centered.hh"

namespace Amanzi {
namespace Operators {

UpwindCellCentered::UpwindCellCentered(std::string pkname,
        std::string cell_coef,
        std::string face_coef) :
    pkname_(pkname),
    cell_coef_(cell_coef),
    face_coef_(face_coef) {};


void UpwindCellCentered::Update(const Teuchos::Ptr<State>& S,
                                  const Teuchos::Ptr<Debugger>& db) {

  const CompositeVector& cell = S->Get<CompositeVector>(cell_coef_);
  CompositeVector& face = S->GetW<CompositeVector>(face_coef_, pkname_);

  *face.ViewComponent("cell") = *cell.ViewComponent("cell");
  if (face.HasComponent("face")) {
    face.ViewComponent("face",true)->PutScalar(1.0);
  }
};


} //namespace
} //namespace
