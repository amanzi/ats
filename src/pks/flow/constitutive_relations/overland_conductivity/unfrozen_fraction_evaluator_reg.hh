/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "unfrozen_fraction_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator,UnfrozenFractionEvaluator>
UnfrozenFractionEvaluator::fac_("unfrozen fraction");

} //namespace
} //namespace

