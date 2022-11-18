
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

  // dummy data
  // 1 col, 100 cells
  double time = 0.;
  int n = 1;
  int m = 100;
  int ncols_local, ncols_global, ncells_per_col;
  std::vector<double> soil_infil(n, 10.0);
  std::vector<double> soil_evap(n, 3.0);
  std::vector<double> air_pres(n);
  std::vector<double> surf_pres(n);
  std::vector<double> elev(n);
  std::vector<double> surf_area_m2(n);
  std::vector<double> lat(n);
  std::vector<double> lon(n);

  std::vector<double> dz(m);
  std::vector<double> depth(m);
  std::vector<double> root_tran(m);
  std::vector<double> soil_pres(m);
  std::vector<double> soil_pot(m);
  std::vector<double> satl(m);
  std::vector<double> sati(m);
  std::vector<int> pft_i(m);

  // dummy fortran comm
  MPI_Fint comm = 0;

  // test driver directly
  auto driver = std::unique_ptr<ATS::ELM_ATSDriver>(ATS::createELM_ATSDriver(&comm, input_filename.data()));
  driver->setup();
  driver->get_mesh_info(ncols_local, ncols_global, lat.data(), lon.data(), elev.data(), surf_area_m2.data(), pft_i.data(), ncells_per_col, depth.data());
  driver->initialize(time, air_pres.data(), soil_pres.data());
  driver->set_potential_sources(soil_infil.data(), soil_evap.data(), root_tran.data());
  driver->advance_test();
  driver->get_waterstate(surf_pres.data(), soil_pres.data(), soil_pot.data(), satl.data(), sati.data());
  driver->finalize();

  // test api
  auto driver_api = ats_create_c(&comm, input_filename.data());
  ats_setup_c(driver_api);
  ats_get_mesh_info_c(driver_api, &ncols_local, &ncols_global, lat.data(), lon.data(),
                      elev.data(), surf_area_m2.data(), pft_i.data(), &ncells_per_col, depth.data());
  ats_initialize_c(driver_api, &time, air_pres.data(), soil_pres.data());
  ats_set_potential_sources_c(driver_api, soil_infil.data(), soil_evap.data(), root_tran.data());
  ats_advance_test_c(driver_api);
  ats_get_waterstate_c(driver_api, surf_pres.data(), soil_pres.data(), soil_pot.data(), satl.data(), sati.data());
  ats_delete_c(driver_api);
  std::cout << "DONE WITH ELM-ATS C++ TEST" << std::endl;

  return 0;
}
