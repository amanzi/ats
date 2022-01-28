/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "carbon_decomposition_rate_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,CarbonDecomposeRateEvaluator> CarbonDecomposeRateEvaluator::reg_("carbon decomposition rate");

}
}
