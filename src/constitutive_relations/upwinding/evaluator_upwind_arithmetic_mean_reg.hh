/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "evaluator_upwind_arithmetic_mean.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<Evaluator, EvaluatorUpwindArithmeticMean>
  EvaluatorUpwindArithmeticMean::factory_("upwind: arithmetic mean");

} // namespace Relations
} // namespace Amanzi
