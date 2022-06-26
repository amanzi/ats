#include "rel_perm_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,RelPermEvaluator> RelPermEvaluator::factory_("WRM rel perm");

} //namespace
} //namespace
