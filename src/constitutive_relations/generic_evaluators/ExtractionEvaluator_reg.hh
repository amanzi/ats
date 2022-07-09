/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
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
Utils::RegisteredFactory<Evaluator,ExtractionEvaluator> ExtractionEvaluator::reg_("extraction evaluator");

} //namespace
} //namespace

