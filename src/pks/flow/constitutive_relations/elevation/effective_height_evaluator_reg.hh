/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluator for determining height( rho, head )

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "effective_height_model.hh"
#include "effective_height_evaluator.hh"


namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,EffectiveHeightEvaluator> EffectiveHeightEvaluator::factory_("effective height");

} //namespace
} //namespace
