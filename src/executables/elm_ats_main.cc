
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
#include "simulation_driver.hh"

// registration files
#include "state_evaluators_registration.hh"

#include "ats_relations_registration.hh"
#include "ats_transport_registration.hh"
#include "ats_energy_pks_registration.hh"
#include "ats_energy_relations_registration.hh"
#include "ats_flow_pks_registration.hh"
#include "ats_flow_relations_registration.hh"
#include "ats_deformation_registration.hh"
#include "ats_bgc_registration.hh"
#include "ats_surface_balance_registration.hh"
#include "ats_mpc_registration.hh"
//#include "ats_sediment_transport_registration.hh"
#include "mdm_transport_registration.hh"
#include "multiscale_transport_registration.hh"
#ifdef ALQUIMIA_ENABLED
#include "pks_chemistry_registration.hh"
#endif


// include fenv if it exists
#include "boost/version.hpp"
#if (BOOST_VERSION / 100 % 1000 >= 46)
#include "boost/config.hpp"
#ifndef BOOST_NO_FENV_H
#ifdef _GNU_SOURCE
#define AMANZI_USE_FENV
#include "boost/detail/fenv.hpp"
#endif
#endif
#endif

#include "boost/filesystem.hpp"

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
  std::cout << "START WITH ELM-ATS C++ TEST" << std::endl;

  // dummy data
  // 1 col, 15 cells
  int n = 1;
  int m = 15;
  //unit: kgH2O/m2/s, + in (sink), - out (source)
  std::vector<double> soil_infil(n, 0.0); // (n, 2.5e-4) - for 'test_interface_infiltration.xml'
  std::vector<double> soil_evap(n,-0.0e-9);
  std::vector<double> root_tran(m,-1.0e-5);   // source: - (out of subsurf), sink: + (into subsurf)

  int ncols_local, ncols_global, ncells_per_col;
  std::vector<double> surf_pres(n);
  std::vector<double> elev(n);
  std::vector<double> surf_area_m2(n);
  std::vector<double> lat(n);
  std::vector<double> lon(n);

  std::vector<double> dz(m);
  std::vector<double> depth(m);
  std::vector<double> soil_pres(m);
  std::vector<double> satl(m);

  // dummy fortran comm
  MPI_Fint comm = 0;

  // test driver directly
  //auto driver = std::make_unique<ATS::ELM_ATSDriver>();
  //driver->setup(&comm, input_filename.data());
  //driver->get_mesh_info(&ncols_local, &ncols_global, &ncells_per_col, dz.data(), depth.data(),
  //  elev.data(), surf_area_m2.data(), lat.data(), lon.data());
  //driver->initialize();
  //driver->set_sources(soil_infil.data(), soil_evap.data(), root_tran.data(), &n, &m);
  //driver->advance_test();
  //driver->get_waterstate(surf_pres.data(), soil_pres.data(), satl.data(), &n, &m);
  //driver->finalize();

  // test api
  auto driver_api = ats_create();
  ats_setup(driver_api, &comm, input_filename.data());
  ats_get_mesh_info(driver_api, &ncols_local, &ncols_global, &ncells_per_col, dz.data(), depth.data(),
    elev.data(), surf_area_m2.data(), lat.data(), lon.data());
  ats_initialize(driver_api);
  ats_set_sources(driver_api, soil_infil.data(), soil_evap.data(), root_tran.data(), &n, &m);
  ats_advance_test(driver_api);
  ats_get_waterstate(driver_api, surf_pres.data(), soil_pres.data(), satl.data(), &n, &m);
  ats_delete(driver_api);
  
  std::cout << "DONE WITH ELM-ATS C++ TEST" << std::endl;

  return 0;
}