/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  An evaluator for pulling the darcy flux, at the surface, from the
  subsurface field and putting it into a surface field.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "overland_source_from_subsurface_flux_evaluator.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<Evaluator, OverlandSourceFromSubsurfaceFluxEvaluator>
  OverlandSourceFromSubsurfaceFluxEvaluator::fac_("overland source from subsurface via flux");

} // namespace Relations
} // namespace Amanzi
