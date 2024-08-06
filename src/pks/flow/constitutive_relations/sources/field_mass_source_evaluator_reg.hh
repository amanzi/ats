/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Phong V.V. Le
*/

#include "field_mass_source_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, FieldMassSourceEvaluator>
  FieldMassSourceEvaluator::reg_("field mass source");

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
