/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  MultiplicativeEvaluator is the generic evaluator for multipying two vectors.

*/

#include "MultiplicativeEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, MultiplicativeEvaluator> MultiplicativeEvaluator::factory_(
  "multiplicative evaluator");

} // namespace Relations
} // namespace Amanzi
