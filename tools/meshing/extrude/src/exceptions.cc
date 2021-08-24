#include "exceptions.hh"

#include <sstream>

namespace Exceptions {

void set_exception_behavior_raise() { behavior = Exceptions::RAISE; }
void set_exception_behavior_abort() { behavior = Exceptions::ABORT; }
void set_exception_behavior(Exception_action action) { behavior = action; }

Exception_action exception_behavior()
{
  return behavior;
}

} // expections error
