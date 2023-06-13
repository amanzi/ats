/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include <UnitTest++.h>

#include <iostream>

#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "AmanziComm.hh"
#include "ats_mesh_factory.hh"

using namespace Amanzi;

struct Runner {
  Runner() {}

  void setup(const std::string& filename)
  {
    std::cout << "Test: " << filename << std::endl;
    comm = getDefaultComm();
    plist = Teuchos::getParametersFromXmlFile(filename);
    gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, plist->sublist("regions"), *comm));
    S = Teuchos::rcp(new State(plist->sublist("state")));
  }

  void go() { ATS::createMeshes(plist->sublist("mesh"), comm, gm, *S); }

  Teuchos::RCP<Teuchos::ParameterList> plist;
  Comm_ptr_type comm;
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<State> S;
};


SUITE(ATS_MESH_FACTORY)
{
  TEST_FIXTURE(Runner, EXTRACT_SURFACE)
  {
    setup("test/executable_mesh_extract_surface.xml");
    go();

    CHECK(S->HasMesh(""));
    CHECK(S->HasMesh("domain"));

    int has_surface = 0;
    int set_size = S->GetMesh("domain")->getSetSize(
      "surface", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    if (set_size > 0) {
      has_surface += 1;
      CHECK(S->HasMesh("surface"));
      CHECK_EQUAL(set_size,
                  S->GetMesh("surface")->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                      AmanziMesh::Parallel_kind::OWNED));
    }

    int total_has_surface = 0;
    comm->SumAll(&has_surface, &total_has_surface, 1);
    CHECK(total_has_surface > 0);
  }
}
