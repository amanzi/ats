/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see ATS_DIR/COPYRIGHT
Author: Ethan Coon (ecoon@lanl.gov)

------------------------------------------------------------------------- */

#include <iostream>

#include <Epetra_MpiComm.h>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "GlobalVerbosity.hh"
#include "VerboseObject.hh"
#include "AmanziComm.hh"
#include "AmanziTypes.hh"

#include "GeometricModel.hh"
#include "coordinator.hh"
#include "State.hh"

#include "errors.hh"
#include "exceptions.hh"

#include "ats_mesh_factory.hh"
#include "simulation_driver.hh"


int SimulationDriver::Run(
    const Amanzi::Comm_ptr_type& comm,
    const Teuchos::RCP<Teuchos::ParameterList>& plist) {

  // verbosity settings
  setDefaultVerbLevel(Amanzi::VerbosityLevel::level_);
  Teuchos::EVerbosityLevel verbLevel = getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = getOStream();
  Teuchos::OSTab tab = getOSTab(); // This sets the line prefix and adds one tab

  // size, rank
  int rank = comm->getRank();
  int size = comm->getSize();

  // print header material
  if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_LOW,true)) {
    // print parameter list
    *out << "======================> dumping parameter list <======================" <<
      std::endl;
    Teuchos::writeParameterListToXmlOStream(*plist, *out);
    *out << "======================> done dumping parameter list. <================" <<
      std::endl;
  }

  // create the geometric model and regions
  auto reg_params = Teuchos::sublist(plist, "regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, *reg_params, *comm) );

  // Create the state.
  auto state_plist = Teuchos::sublist(plist, "state");
  auto S = Teuchos::rcp(new Amanzi::State(state_plist));

  // create and register meshes
  ATS::createMeshes(*plist, comm, gm, *S);
 
  // create the top level Coordinator
  ATS::Coordinator coordinator(plist, S, comm);
  
  // run the simulation
  coordinator.CycleDriver();
  return 0;
}


