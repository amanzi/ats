/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "ComponentAssemblerEvaluator.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<Evaluator, ComponentAssemblerEvaluator>
  ComponentAssemblerEvaluator::factory_("component assembler");

} // namespace Relations
} // namespace Amanzi
