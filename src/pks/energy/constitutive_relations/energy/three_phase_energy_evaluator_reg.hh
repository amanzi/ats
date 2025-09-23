/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "three_phase_energy_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {
namespace Relations {

Utils::RegisteredFactory<Evaluator, ThreePhaseEnergyEvaluator> ThreePhaseEnergyEvaluator::reg_(
  "three phase energy");

} // namespace Relations
} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
