#include "eos_ideal_gas_evaluator.hh"

namespace Amanzi {
namespace General {
namespace Relations {

Utils::RegisteredFactory<Evaluator,EosIdealGasEvaluator> EosIdealGasEvaluator::reg_("ideal gas equation of state");

} //namespace
} //namespace
} //namespace
