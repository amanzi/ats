/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "volumetric_snow_ponded_depth_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, VolumetricSnowPondedDepthEvaluator>
  VolumetricSnowPondedDepthEvaluator::reg_("volumetric ponded and snow depths");

} // namespace Flow
} // namespace Amanzi
