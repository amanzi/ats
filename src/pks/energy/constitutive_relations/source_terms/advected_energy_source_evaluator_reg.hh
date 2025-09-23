/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Source term evaluator for enthalpy of mass source.

*/

#include "advected_energy_source_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

Utils::RegisteredFactory<Evaluator, AdvectedEnergySourceEvaluator>
  AdvectedEnergySourceEvaluator::factory_("advected energy source");

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
