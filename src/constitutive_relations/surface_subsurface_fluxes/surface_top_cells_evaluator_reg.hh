/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Specifies a value on the surface from the value in the cell just below the
  surface.

*/

#include "surface_top_cells_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, SurfaceTopCellsEvaluator> SurfaceTopCellsEvaluator::reg_(
  "surface from top cell evaluator");

} // namespace Relations
} // namespace Amanzi
