/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  A field evaluator for an unchanging cell volume.
*/

#include "EvaluatorSubgridReturn.hh"

namespace Amanzi {

Utils::RegisteredFactory<Evaluator, EvaluatorSubgridReturn> EvaluatorSubgridReturn::fac_(
  "subgrid return");

} // namespace Amanzi
