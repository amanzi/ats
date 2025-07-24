/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include "dbc.hh"

#include <sstream>

namespace DBC {

Assertion::Assertion(const char* assertion, const char* filename, unsigned int line_number)
  : assertion_(assertion), filename_(filename), line_number_(line_number) {};


const char*
Assertion::what() const throw()
{
  std::ostringstream message;
  message << "Assertion: \"" << assertion_ << "\" failed in file: " << filename_
          << ", at line: " << line_number_ << std::endl;

  return message.str().c_str();
}


void
amanzi_assert(const char* cond, const char* file, unsigned int line)
{
  Exceptions::amanzi_throw(Assertion(cond, file, line));
}

} // namespace DBC
