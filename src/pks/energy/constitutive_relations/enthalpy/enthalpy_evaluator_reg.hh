/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "enthalpy_evaluator.hh"

namespace Amanzi {
namespace Energy {

// registry of method
Utils::RegisteredFactory<Evaluator, EnthalpyEvaluator> EnthalpyEvaluator::factory_("enthalpy");

} // namespace Energy
} // namespace Amanzi
