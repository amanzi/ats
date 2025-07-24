/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "eos_ideal_gas_evaluator.hh"

namespace Amanzi {
namespace General {
namespace Relations {

Utils::RegisteredFactory<Evaluator, EosIdealGasEvaluator> EosIdealGasEvaluator::reg_(
  "ideal gas equation of state");

} // namespace Relations
} // namespace General
} // namespace Amanzi
