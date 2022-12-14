/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "exceptions.hh"

#include <sstream>

namespace Exceptions {

void
set_exception_behavior_raise()
{
  behavior = Exceptions::RAISE;
}
void
set_exception_behavior_abort()
{
  behavior = Exceptions::ABORT;
}
void
set_exception_behavior(Exception_action action)
{
  behavior = action;
}

Exception_action
exception_behavior()
{
  return behavior;
}

} // namespace Exceptions
