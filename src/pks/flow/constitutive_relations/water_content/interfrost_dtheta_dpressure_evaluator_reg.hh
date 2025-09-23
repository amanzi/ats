/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "interfrost_dtheta_dpressure_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, InterfrostDthetaDpressureEvaluator>
  InterfrostDthetaDpressureEvaluator::reg_("interfrost dtheta_dpressure");

} // namespace Relations
} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
