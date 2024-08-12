#include <cstdlib>
#include <iostream>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_simple.hh"
#include "State.hh"

#include "lake_thermo_pk.hh"

#include "lake_thermo_test_class.hh"

/* **************************************************************** */
TEST(ADVANCE_WITH_SIMPLE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: advance using simple mesh" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  // read parameter list
  Teuchos::ParameterList parameter_list;
  string xmlFileName = "test/two_phase_advance_simple.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework
  Teuchos::ParameterList region_list =
    parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
  Teuchos::RCP<Mesh> mesh = Teuchos::rcp(new Mesh_simple(0.0,0.0,0.0, 1.0,1.0,1.0, 20, 20, 2, comm, gm));

  // create and initialize the test class
  EnergyTestOne test(parameter_list, mesh, 1);
  test.initialize();

  // advance the state
  int iter, k;
  double T = 0.0;

  Teuchos::RCP<const CompositeVector> temp =
    test.S1->GetFieldData("temperature");

  iter = 0;
  if (iter < 10) {
    printf( "T=%6.2f  C_0(x):", T );
    for( int k=0; k<15; k++ ) printf("%7.4f", (*temp)(0,k));
    cout << endl;
  }

  double L1, L2;
  while (T < 1.0) {
    double dT = test.EPK->get_dt();
    test.EPK->advance(dT);
    T += dT;
    iter++;

    if (iter < 10) {
      printf( "T=%6.2f  C_0(x):", T );
      for( int k=0; k<15; k++ ) printf("%7.4f", (*temp)(0,k));
      cout << endl;
    }

    test.evaluate_error_temp(T, &L1, &L2);
    CHECK_CLOSE(0., L2, 1e-6);
    CHECK_CLOSE(0., L1, 1e-6);
    test.commit_step();
  }
  delete comm;
}
