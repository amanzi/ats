
#include <iostream>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "VerboseObject_objs.hh"
#include "VerboseObject.hh"

#include "ats_version.hh"
#include "amanzi_version.hh"
#include "tpl_versions.h"

#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "GeometricModel.hh"
#include "State.hh"

#include "errors.hh"
#include "exceptions.hh"
#include "dbc.hh"

#include "ats_mesh_factory.hh"
#include "elm_ats_coordinator.hh"
#include "elm_ats_driver.hh"

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

namespace ATS {

int
ELM_ATSDriver::setup(char *infile)
{
  // -- create communicator & get process rank
  auto comm = Amanzi::getDefaultComm();
  auto rank = comm->MyPID();

  // convert input file to std::string for easier handling
  // infile must be null-terminated
  std::string input_filename(infile);

  // parse the input file and check validity
  if (input_filename.empty()) {
    if (rank == 0) {
      std::cerr << "ERROR: no input file provided" << std::endl;
    }
    return 1;
  } else if (!boost::filesystem::exists(input_filename)) {
    if (rank == 0) {
      std::cerr << "ERROR: input file \"" << input_filename << "\" does not exist." << std::endl;
    }
    return 1;
  }

  // -- parse input file
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(input_filename);

  // -- set default verbosity level to no output
  Amanzi::VerboseObject::global_default_level = Teuchos::VERB_NONE;

  // create the geometric model and regions
  Teuchos::ParameterList reg_params = plist->sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, *comm) );

  // Create the state.
  Teuchos::ParameterList state_plist = plist->sublist("state");
  S_ = Teuchos::rcp(new Amanzi::State(state_plist));

  // create and register meshes
  ATS::Mesh::createMeshes(*plist, comm, gm, *S_);

  // create ELM coordinator object
  elm_coordinator = std::make_unique<ELM_ATSCoordinator>(*plist, S_, comm);
  // call coordinator setup
  elm_coordinator->setup();

  return 0;
}

void ELM_ATSDriver::initialize()
{
  elm_coordinator->initialize();
}


void ELM_ATSDriver::advance(double *dt)
{
  elm_coordinator->advance(*dt);
}


} // namespace ATS
