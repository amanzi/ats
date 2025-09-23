/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/* -------------------------------------------------------------------------

   ATS

   Two-phase TC evaluator
   ------------------------------------------------------------------------- */

#include "thermal_conductivity_twophase_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

// registry of method
Utils::RegisteredFactory<Evaluator, ThermalConductivityTwoPhaseEvaluator>
  ThermalConductivityTwoPhaseEvaluator::factory_("two-phase thermal conductivity");

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
