/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Extracts a field on one mesh from a field on a superset of that mesh using
parent entities.

*/

#include "ExtractionEvaluator.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<Evaluator, ExtractionEvaluator> ExtractionEvaluator::reg_(
  "extraction evaluator");

} // namespace Relations
} // namespace Amanzi
