/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

#include "evaporation_downregulation_lee_pielke_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, EvaporationDownregulationLeePielkeEvaluator>
  EvaporationDownregulationLeePielkeEvaluator::fac_("evaporation downregulation, Lee-Pielke");

} // namespace Flow
} // namespace Amanzi
