#include "wrm_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,WRMEvaluator> WRMEvaluator::factory_("WRM");
Utils::RegisteredFactory<Evaluator,WRMEvaluator> WRMEvaluator::factory2_("wrm");

} //namespace
} //namespace
