#include "enthalpy_evaluator.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<Evaluator, EnthalpyEvaluator> EnthalpyEvaluator::factory_("enthalpy");

} // namespace Energy
} // namespace Amanzi
