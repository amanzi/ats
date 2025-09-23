/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "liquid_gas_water_content_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, LiquidGasWaterContentEvaluator>
  LiquidGasWaterContentEvaluator::reg_("liquid+gas water content");

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
