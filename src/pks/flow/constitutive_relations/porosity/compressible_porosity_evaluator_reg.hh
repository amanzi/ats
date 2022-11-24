/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the porosity, given a small compressibility of rock.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "compressible_porosity_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, CompressiblePorosityEvaluator>
  CompressiblePorosityEvaluator::fac_("compressible porosity");

} // namespace Flow
} // namespace Amanzi
