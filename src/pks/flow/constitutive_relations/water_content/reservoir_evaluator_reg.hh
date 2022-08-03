#include "reservoir_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator,ReservoirEvaluator> ReservoirEvaluator::reg_("reservoir model");

} //namespace
} //namespace
} //namespace
