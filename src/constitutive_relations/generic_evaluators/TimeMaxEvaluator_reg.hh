/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  TimeMaxEvaluator is the generic evaluator for multipying two vectors.

*/

#include "TimeMaxEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, TimeMaxEvaluator> TimeMaxEvaluator::reg_("max in time");

} // namespace Relations
} // namespace Amanzi
