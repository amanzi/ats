/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "surface_distributed_tiles_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, SurfDistributedTilesRateEvaluator>
  SurfDistributedTilesRateEvaluator::reg_("surface distributed tiles");

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
