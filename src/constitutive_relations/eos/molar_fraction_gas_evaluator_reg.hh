/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "vapor_pressure_relation_factory.hh"
#include "molar_fraction_gas_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator,MolarFractionGasEvaluator> MolarFractionGasEvaluator::factory_("molar fraction gas");

} // namespace
} // namespace

