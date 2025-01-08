/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Specifies a value on the surface from the value in the cell just below the
  surface.

*/

#include "top_cells_surface_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, TopCellsSurfaceEvaluator>
  TopCellsSurfaceEvaluator::reg_("top cell from surface evaluator");

} // namespace Relations
} // namespace Amanzi
