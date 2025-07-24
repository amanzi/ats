/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "elevation_evaluator_column.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, ColumnElevationEvaluator> ColumnElevationEvaluator::reg_(
  "column elevation");

} // namespace Flow
} // namespace Amanzi
