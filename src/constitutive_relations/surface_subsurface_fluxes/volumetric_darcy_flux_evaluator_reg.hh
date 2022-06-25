
/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  An evaluator for converting the darcy flux to volumetric flux

  Authors: Daniil Svyatsky  (dasvyat@lanl.gov)
*/


#include "volumetric_darcy_flux_evaluator.hh"

namespace Amanzi {
namespace Relations {

Utils::RegisteredFactory<Evaluator,Volumetric_FluxEvaluator>
Volumetric_FluxEvaluator::fac_("volumetric darcy flux");

} // namespace
} // namespace
