/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "surface_culvert_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, SurfCulvertEvaluator> SurfCulvertEvaluator::reg_("culvert");

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
