/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ATS version) (coonet@ornl.gov)
*/

/* -------------------------------------------------------------------------
  A high level test class for the MatrixMFD operator.

------------------------------------------------------------------------- */

#include "divgrad_test.hh"

namespace Amanzi {
namespace TestPKs {

RegisteredPKFactory_ATS<DivGradTest> DivGradTest::reg_("div-grad operator test");

} // namespace TestPKs
} // namespace Amanzi
