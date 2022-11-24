/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
  A high level test class for the MatrixMFD operator.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */

#include "divgrad_test.hh"

namespace Amanzi {
namespace TestPKs {

RegisteredPKFactory_ATS<DivGradTest> DivGradTest::reg_("div-grad operator test");

} // namespace TestPKs
} // namespace Amanzi
