/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "column_average_temp_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,ColumnAverageTempEvaluator> ColumnAverageTempEvaluator::reg_("column average temperature");

}
}
