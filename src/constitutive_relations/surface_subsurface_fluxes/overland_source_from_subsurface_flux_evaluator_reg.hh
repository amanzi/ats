/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

*/

#include "overland_source_from_subsurface_flux_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

Utils::RegisteredFactory<Evaluator, OverlandSourceFromSubsurfaceFluxEvaluator>
  OverlandSourceFromSubsurfaceFluxEvaluator::fac_("overland source from subsurface via flux");

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
