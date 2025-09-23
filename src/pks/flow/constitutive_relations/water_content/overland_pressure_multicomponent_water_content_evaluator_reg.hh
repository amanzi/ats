/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "overland_pressure_multicomponent_water_content_evaluator.hh"
namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
Utils::RegisteredFactory<Evaluator, OverlandPressureMulticomponentWaterContentEvaluator>
  OverlandPressureMulticomponentWaterContentEvaluator::reg_(
    "overland pressure multicomponent water content");
}
} // namespace ATS_Physics
} // namespace Amanzi
