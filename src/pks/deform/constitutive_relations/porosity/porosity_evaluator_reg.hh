/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the porosity after a cell volume change due to mesh deformation.

  Authors: Markus Berndt (berndt@lanl.gov)
*/

#include "porosity_evaluator.hh"

namespace Amanzi {
namespace Deform {
namespace DeformRelations {

Utils::RegisteredFactory<Evaluator,PorosityEvaluator> PorosityEvaluator::factory_("deformation porosity");

} //namespace
} //namespace
} //namespace

