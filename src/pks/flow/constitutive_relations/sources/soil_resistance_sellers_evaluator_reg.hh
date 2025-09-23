/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Bo Gao (gaob@ornl.gov)
*/

#include "soil_resistance_sellers_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

// registry of method
Utils::RegisteredFactory<Evaluator, SoilResistanceSellersEvaluator>
  SoilResistanceSellersEvaluator::fac_("soil resistance, Sellers");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
