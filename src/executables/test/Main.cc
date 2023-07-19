/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <mpi.h>

#include <TestReporterStdout.h>
#include "Teuchos_GlobalMPISession.hpp"
#include <UnitTest++.h>

#include "ats_registration_files.hh"
#include "VerboseObject_objs.hh"

#include "Kokkos_Core.hpp"

int
main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv); 
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  auto result = UnitTest::RunAllTests();
  Kokkos::finalize();
  return result; 
}
