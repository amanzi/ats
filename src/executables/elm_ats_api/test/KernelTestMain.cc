/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rich Fiorella, Ethan Coon
*/

//! Minimal UnitTest++ runner for the dependency-free water-balance kernel test.
/*!
  The kernel (elm_ats_water_balance_error_c) uses only iso_c_binding, so this
  runner deliberately avoids the Kokkos/MPI initialization done in the shared
  executables/test/Main.cc -- none of it is needed here.
*/

#include <UnitTest++.h>

int
main(int /*argc*/, char* /*argv*/[])
{
  return UnitTest::RunAllTests();
}
