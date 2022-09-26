/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Specifies a value on the surface from the value in the cell just below the
  surface.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "surface_top_cells_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,SurfaceTopCellsEvaluator> SurfaceTopCellsEvaluator::reg_("surface from top cell evaluator");

} //namespace
} //namespace

