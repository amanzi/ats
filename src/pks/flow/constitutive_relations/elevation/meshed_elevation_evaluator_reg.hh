/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,MeshedElevationEvaluator> MeshedElevationEvaluator::reg_("meshed elevation");
Utils::RegisteredFactory<Evaluator,StandaloneElevationEvaluator> StandaloneElevationEvaluator::reg_("standalone elevation");

}
}
