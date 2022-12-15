/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/

#include "carbon_decomposition_rate_evaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, CarbonDecomposeRateEvaluator>
  CarbonDecomposeRateEvaluator::reg_("carbon decomposition rate");

} // namespace Relations
} // namespace Amanzi
