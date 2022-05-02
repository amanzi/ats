/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see ATS_DIR/COPYRIGHT
Author: Ethan Coon (ecoon@lanl.gov)

------------------------------------------------------------------------- */

#include <iostream>

#include <Epetra_MpiComm.h>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"

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

namespace ATS {

int
SimulationDriver::Run(const Teuchos::RCP<const Amanzi::Comm_type>& comm,
                      Teuchos::ParameterList& plist)
{
  Amanzi::VerboseObject vo("Simulation Driver", plist);
  Teuchos::OSTab tab = vo.getOSTab();

  // print header material
  if (vo.os_OK(Teuchos::VERB_LOW)) {
    // print parameter list
    *vo.os() << "======================> dumping parameter list <======================" <<
      std::endl;
    Teuchos::writeParameterListToXmlOStream(plist, *vo.os());
    *vo.os() << "======================> done dumping parameter list. <================" <<
      std::endl;
  }

  // create the geometric model and regions
  Teuchos::ParameterList reg_params = plist.sublist("regions");
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
    Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, *comm) );

  // Create the state.
  Teuchos::ParameterList state_plist = plist.sublist("state");
  Teuchos::RCP<Amanzi::State> S = Teuchos::rcp(new Amanzi::State(state_plist));

  // create and register meshes
  //ATS::createMeshes(plist.sublist("mesh"), comm, gm, *S);
  ATS::Mesh::createMeshes(plist, comm, gm, *S);

  // create the top level Coordinator
  ATS::Coordinator coordinator(plist, S, comm);

  // run the simulation
  coordinator.cycle_driver();
  return 0;
}

} // namespace
