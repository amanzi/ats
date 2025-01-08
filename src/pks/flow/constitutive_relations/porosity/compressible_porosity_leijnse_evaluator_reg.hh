/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Evaluates the porosity, given a small compressibility of rock.

*/

#include "compressible_porosity_leijnse_evaluator.hh"

namespace Amanzi {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, CompressiblePorosityLeijnseEvaluator>
  CompressiblePorosityLeijnseEvaluator::fac_("compressible porosity leijnse");

} // namespace Flow
} // namespace Amanzi
