#include "thaw_depth_columns_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,ThawDepthColumnsEvaluator> ThawDepthColumnsEvaluator::reg_("thaw depth, columns");

}
}
