/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ViscosityEvaluator is the interface between state/data and the model, a VPM.

*/

#include "viscosity_relation_factory.hh"
#include "viscosity_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, ViscosityEvaluator> ViscosityEvaluator::factory_("viscosity");

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
