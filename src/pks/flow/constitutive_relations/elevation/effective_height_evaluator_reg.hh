/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "effective_height_model.hh"
#include "effective_height_evaluator.hh"


namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, EffectiveHeightEvaluator>
  EffectiveHeightEvaluator::factory_("effective height");

} // namespace Flow
} // namespace Amanzi
