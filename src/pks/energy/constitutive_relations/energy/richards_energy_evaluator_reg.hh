/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "richards_energy_evaluator.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

Utils::RegisteredFactory<Evaluator, RichardsEnergyEvaluator> RichardsEnergyEvaluator::reg_(
  "richards energy");

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
