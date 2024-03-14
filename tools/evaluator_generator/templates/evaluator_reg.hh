/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "{evalName}_evaluator.hh"

namespace Amanzi {
{
  namespace {namespace}
  {
    {
      namespace Relations {
      {
        Utils::RegisteredFactory<Evaluator, { evalClassName } Evaluator>{
          evalClassName
        } Evaluator::reg_("{evalNameString}");

      }
      } // namespace Relations
    }
  } // namespace }
} // namespace Amanzi
} // namespace Amanzi
