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

#include "dbc.hh"
#include "errors.hh"
#include "simulation_driver.hh"

#include "state_evaluators_registration.hh"

#include "constitutive_relations_eos_registration.hh"
#include "constitutive_relations_surface_subsurface_fluxes_registration.hh"
#include "constitutive_relations_generic_evaluators_registration.hh"

#include "flow_relations_registration.hh"
#include "flow_icy_overland_registration.hh"
#include "flow_overland_pressure_registration.hh"
#include "flow_permafrost_registration.hh"
#include "flow_richards_registration.hh"
#include "transport_amanzi_registration.hh"
//#include "pks_transport_registration.hh"
#include "multiscale_transport_registration.hh"
#include "mdm_transport_registration.hh"

#include "surface_balance_SEB_registration.hh"
#include "BGC_registration.hh"

//#include "chemistry_amanzi_registration.hh"
//#include "pks_chemistry_registration.hh"
#include "mpc_pks_registration.hh"
#include "mpc_registration.hh"

#include "deformation_registration.hh"

#include "energy_advection_diffusion_registration.hh"
// #include "energy_constant_temperature_registration.hh"
#include "energy_relations_registration.hh"
#include "energy_surface_ice_registration.hh"
#include "energy_two_phase_registration.hh"
#include "energy_three_phase_registration.hh"






// include fenv if it exists
#include "boost/version.hpp"
#if (BOOST_VERSION / 100 % 1000 >= 46)
 #include "boost/config.hpp"
 #ifndef BOOST_NO_FENV_H
//  #ifdef _GNU_SOURCE
   #define AMANZI_USE_FENV
   #include "boost/detail/fenv.hpp"
//  #endif
 #endif
#endif


Teuchos::EVerbosityLevel Amanzi::VerbosityLevel::level_ = Teuchos::VERB_MEDIUM;

int main(int argc, char *argv[])
{

#ifdef AMANZI_USE_FENV
  //  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  feraiseexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif

  Teuchos::GlobalMPISession mpiSession(&argc,&argv,0);

  Teuchos::CommandLineProcessor CLP;
  CLP.setDocString("\nATS: simulations for ecosystem hydrology\n");

  std::string xmlInFileName = "options.xml";
  CLP.setOption("xml_file", &xmlInFileName, "XML options file");
  CLP.throwExceptions(false);
  
  Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parseReturn = CLP.parse(argc, argv);

  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    return 0;
  }
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return 1;
  }

  MPI_Comm mpi_comm(MPI_COMM_WORLD);

  // read the main parameter list
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xmlInFileName); 

  Teuchos::RCP<Teuchos::FancyOStream> fos;
  Teuchos::readVerboseObjectSublist(&*plist, &fos, &Amanzi::VerbosityLevel::level_);

  SimulationDriver simulator;
  int ret = simulator.Run(mpi_comm, *plist);
}


