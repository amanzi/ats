/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Phong V.V. Le (lepv@ornl.gov)
*/

#include "qc_relation_multi_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

Utils::RegisteredFactory<Evaluator, QCRelationMultiEvaluator>
  QCRelationMultiEvaluator::reg_("q-c multi");

} // namespace Relations
} // namespace Flow
} // namespace Amanzi
