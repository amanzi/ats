/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky  (dasvyat@lanl.gov)
*/

/*
  An evaluator for converting the darcy flux to volumetric flux

*/


#include "volumetric_darcy_flux_evaluator.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<Evaluator, Volumetric_FluxEvaluator> Volumetric_FluxEvaluator::fac_(
  "volumetric darcy flux");

} // namespace Relations
} // namespace Amanzi
