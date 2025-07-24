/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ViscosityEvaluator is the interface between state/data and the model, a VPM.

*/

#include "pc_ice_water.hh"
#include "pc_ice_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, PCIceEvaluator> PCIceEvaluator::factory_(
  "capillary pressure, water over ice");

} // namespace Flow
} // namespace Amanzi
