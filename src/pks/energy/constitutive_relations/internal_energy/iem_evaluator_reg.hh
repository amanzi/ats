/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

*/

#include "iem_evaluator.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<Evaluator, IEMEvaluator> IEMEvaluator::factory_("iem");


} // namespace Energy
} // namespace Amanzi
