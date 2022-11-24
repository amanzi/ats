#include "ColumnSumEvaluator.hh"
#include "activelayer_average_temp_evaluator.hh"
#include "thaw_depth_evaluator.hh"
#include "water_table_depth_evaluator.hh"


namespace Amanzi {
namespace Relations {

// registry of method
template <>
Utils::RegisteredFactory<Evaluator, ColumnSumEvaluator>
  ColumnSumEvaluator::reg_("column sum evaluator");
template <>
Utils::RegisteredFactory<Evaluator, ActiveLayerAverageTempEvaluator>
  ActiveLayerAverageTempEvaluator::reg_("active layer average temperature");
template <>
Utils::RegisteredFactory<Evaluator, ThawDepthEvaluator> ThawDepthEvaluator::reg_("thaw depth");
template <>
Utils::RegisteredFactory<Evaluator, WaterTableDepthEvaluator>
  WaterTableDepthEvaluator::reg_("water table depth");

} // namespace Relations
} // namespace Amanzi
