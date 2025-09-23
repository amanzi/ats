/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "micropore_macropore_flux_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace SurfaceBalance {
namespace Relations {

Utils::RegisteredFactory<Evaluator, MicroporeMacroporeFluxEvaluator>
  MicroporeMacroporeFluxEvaluator::reg_("micropore-macropore flux");

} // namespace Relations
} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
