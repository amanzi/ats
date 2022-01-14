#include "water_table_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,WaterTableDepthEvaluator> WaterTableDepthEvaluator::reg_("water table depth");

}
}
