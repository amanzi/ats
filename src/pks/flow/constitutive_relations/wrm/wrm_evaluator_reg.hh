#include "wrm_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,WRMEvaluator> WRMEvaluator::factory_("WRM");
Utils::RegisteredFactory<FieldEvaluator,WRMEvaluator> WRMEvaluator::factory2_("wrm");

} //namespace
} //namespace
