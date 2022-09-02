/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Specifies a value on the surface from the value in the cell just below the
  surface.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "top_cells_surface_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,TopCellsSurfaceEvaluator> TopCellsSurfaceEvaluator::reg_("top cell from surface evaluator");

} //namespace
} //namespace

