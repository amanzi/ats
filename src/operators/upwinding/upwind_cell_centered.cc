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

UpwindCellCentered::UpwindCellCentered(const std::string& pkname, const Tag& tag)
  : pkname_(pkname), tag_(tag)
{}


void
UpwindCellCentered::Update(const CompositeVector& cells,
                           CompositeVector& faces,
                           const State& S,
                           const Teuchos::Ptr<Debugger>& db) const
{
  *faces.ViewComponent("cell") = *cells.ViewComponent("cell");
  if (faces.HasComponent("face")) { faces.ViewComponent("face", true)->PutScalar(1.0); }
};


} // namespace Operators
} // namespace Amanzi
