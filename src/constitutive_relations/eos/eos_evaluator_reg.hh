/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  EOSEvaluator is the interface between state/data and the model, an EOS.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "eos_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,EOSEvaluator> EOSEvaluator::factory_("eos");

}
}
