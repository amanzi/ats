#include "pipe_drain_evaluator.hh"
namespace Amanzi {
namespace Flow {
Utils::RegisteredFactory<Evaluator, PipeDrainEvaluator> PipeDrainEvaluator::reg_("pipe drain");
}
} // namespace Amanzi
