/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "microtopographic_relief_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,MicrotopographicReliefEvaluator> MicrotopographicReliefEvaluator::reg_("microtopographic relief");

}
}
