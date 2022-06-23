/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ViscosityEvaluator is the interface between state/data and the model, a VPM.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "pc_ice_water.hh"
#include "pc_ice_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,PCIceEvaluator> PCIceEvaluator::factory_("capillary pressure, water over ice");

} // namespace
} // namespace
