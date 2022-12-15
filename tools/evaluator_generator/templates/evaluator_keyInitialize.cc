/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

// dependency: {arg}
{ var } _key_ = Keys::readKey(plist_, domain_name, "{argString}", "{arg}");
dependencies_.insert({ var } _key_);
