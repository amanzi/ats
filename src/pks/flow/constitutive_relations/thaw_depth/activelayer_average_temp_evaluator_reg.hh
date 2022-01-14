/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "activelayer_average_temp_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,ActiveLayerAverageTempEvaluator> ActiveLayerAverageTempEvaluator::reg_("active layer average temperature");

}
}
