#include "fractional_conductance_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// registry of method
Utils::RegisteredFactory<Evaluator,FractionalConductanceEvaluator> FractionalConductanceEvaluator::factory_("fractional conductance");

}
}
}
