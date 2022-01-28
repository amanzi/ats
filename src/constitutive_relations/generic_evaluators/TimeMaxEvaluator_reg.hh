/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  TimeMaxEvaluator is the generic evaluator for multipying two vectors.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "TimeMaxEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,TimeMaxEvaluator> TimeMaxEvaluator::reg_("max in time");

} // namespace
} // namespace
