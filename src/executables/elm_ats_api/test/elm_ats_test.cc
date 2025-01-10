
#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "VerboseObject_objs.hh"

#include "ats_version.hh"
#include "amanzi_version.hh"
#include "tpl_versions.h"

#include "dbc.hh"
#include "errors.hh"

// registration files
#include "ats_registration_files.hh"


#include "elm_ats_driver.hh"
#include "elm_ats_api.h"

int main(int argc, char *argv[])
{

#ifdef AMANZI_USE_FENV
  //  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  feraiseexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  Teuchos::GlobalMPISession mpiSession(&argc,&argv,0);
  int rank = mpiSession.getRank();

  std::string input_filename;
  if ((argc >= 2) && (argv[argc-1][0] != '-')) {
    input_filename = std::string(argv[argc-1]);
    argc--;
  }

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("Run ATS simulations for ecosystem hydrology.\n\nStandard usage: ats input.xml\n");

  std::string opt_input_filename = "";
  clp.setOption("xml_file", &opt_input_filename, "XML input file");

  bool version(false);
  clp.setOption("version", "no_version", &version, "Print version number and exit.");

  bool print_version(false);
  clp.setOption("print_version", "no_print_version", &print_version, "Print full version info and exit.");

  std::string verbosity;
  clp.setOption("verbosity", &verbosity, "Default verbosity level: \"none\", \"low\", \"medium\", \"high\", \"extreme\".");

  clp.throwExceptions(false);
  clp.recogniseAllOptions(true);

  auto parseReturn = clp.parse(argc, argv);
  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    return 0;
  }
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return 1;
  }

  if (input_filename.empty() && !opt_input_filename.empty()) input_filename = opt_input_filename;

  // dummy data
  double time = 0.;
  int ncols = 5;
  int ncells = 15;
  int ncols_local, ncols_global, ncells_per_col;
  std::vector<double> soil_infil(ncols, 10.0);
  std::vector<double> soil_evap(ncols, 3.0);
  std::vector<double> root_tran(ncols, 6.0);
  std::vector<double> soil_pres(ncols*ncells);
  std::vector<double> satl(ncols*ncells);

  // dummy fortran comm
  MPI_Fint comm = 0;

  // test driver directly
 auto driver = std::unique_ptr<ATS::ELM_ATSDriver>(ATS::createELM_ATSDriver(&comm, input_filename.data()));
 driver->setup();
 driver->initialize(time, soil_pres.data(), satl.data());
 driver->advance_test();
 driver->finalize();
 std::cout << "DONE WITH DRIVER TEST" << std::endl;

  // test api

  auto driver_api = ats_create(&comm, input_filename.data());
  ats_setup(driver_api);
  ats_initialize(driver_api, &time, soil_pres.data(), satl.data());
  ats_advance_test(driver_api);
  ats_delete(driver_api);
  std::cout << "DONE WITH API TEST" << std::endl;

  return 0;
}
