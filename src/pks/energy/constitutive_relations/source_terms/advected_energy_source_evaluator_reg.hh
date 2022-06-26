/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Source term evaluator for enthalpy of mass source.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "advected_energy_source_evaluator.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<Evaluator,AdvectedEnergySourceEvaluator> AdvectedEnergySourceEvaluator::factory_("advected energy source");

} //namespace
} //namespace

