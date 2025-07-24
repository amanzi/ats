/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  AdditiveEvaluator is the generic evaluator for multipying two vectors.

*/

#include "AdditiveEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, AdditiveEvaluator> AdditiveEvaluator::factory_(
  "additive evaluator");

} // namespace Relations
} // namespace Amanzi
