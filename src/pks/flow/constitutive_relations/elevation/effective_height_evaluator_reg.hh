/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluator for determining height( rho, head )

*/

#include "effective_height_model.hh"
#include "effective_height_evaluator.hh"


namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, EffectiveHeightEvaluator> EffectiveHeightEvaluator::factory_(
  "effective height");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
