#include <mpi.h>

#include <TestReporterStdout.h>
#include "Teuchos_GlobalMPISession.hpp"
#include <UnitTest++.h>

#include "ats_registration_files.hh"
#include "VerboseObject_objs.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  return UnitTest::RunAllTests();
}

