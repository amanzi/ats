/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "iem_evaluator.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<Evaluator, IEMEvaluator> IEMEvaluator::factory_("iem");


} // namespace Energy
} // namespace Amanzi
