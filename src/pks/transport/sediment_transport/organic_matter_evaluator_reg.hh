/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  The erosion evaluator gets the erosion rates.

*/

#include "organic_matter_evaluator.hh"

namespace Amanzi {

// registry of method
Utils::RegisteredFactory<Evaluator, OrganicMatterRateEvaluator>
  OrganicMatterRateEvaluator ::factory_("organic matter rate");

} // namespace Amanzi
