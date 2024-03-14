/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, MeshedElevationEvaluator>
  MeshedElevationEvaluator::reg_("meshed elevation");
Utils::RegisteredFactory<Evaluator, StandaloneElevationEvaluator>
  StandaloneElevationEvaluator::reg_("standalone elevation");

} // namespace Flow
} // namespace Amanzi
