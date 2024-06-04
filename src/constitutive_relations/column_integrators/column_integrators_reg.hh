/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "ColumnSumEvaluator.hh"
//#include "activelayer_average_temp_evaluator.hh"
//#include "thaw_depth_evaluator.hh"
#include "water_table_depth_evaluator.hh"


namespace Amanzi {
namespace Relations {

// registry of method
template <>
const std::string ColumnSumEvaluator::eval_type = "column sum evaluator";
template <>
REGISTER(ColumnSumEvaluator);

template <>
const std::string WaterTableDepthEvaluator::eval_type = "water table depth";
template <>
REGISTER(WaterTableDepthEvaluator);

// template <>
// Utils::RegisteredFactory<Evaluator, ActiveLayerAverageTempEvaluator>
//   ActiveLayerAverageTempEvaluator::reg_("active layer average temperature");
// template <>
// Utils::RegisteredFactory<Evaluator, ThawDepthEvaluator> ThawDepthEvaluator::reg_("thaw depth");

} // namespace Relations
} // namespace Amanzi
