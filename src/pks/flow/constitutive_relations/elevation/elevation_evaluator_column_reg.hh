/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "elevation_evaluator_column.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,ColumnElevationEvaluator> ColumnElevationEvaluator::reg_("column elevation");

}
}
