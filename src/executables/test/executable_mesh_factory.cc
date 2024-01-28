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
    std::cout << std::endl
              << std::endl
              << "Test: " << filename << std::endl
              << "---------------------------------------" << std::endl;
    comm = getDefaultComm();
    plist = Teuchos::getParametersFromXmlFile(filename);
    gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, plist->sublist("regions"), *comm));
    S = Teuchos::rcp(new State(plist->sublist("state")));
  }

  void go() { ATS::Mesh::createMeshes(*plist, comm, gm, *S); }

  Teuchos::RCP<Teuchos::ParameterList> plist;
  Comm_ptr_type comm;
  Teuchos::RCP<AmanziGeometry::GeometricModel> gm;
  Teuchos::RCP<State> S;
};


SUITE(ATS_MESH_FACTORY)
{
#if 0
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
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &has_surface, &total_has_surface);
    CHECK(total_has_surface > 0);
  }
#endif

  TEST_FIXTURE(Runner, EXTRACT_SUBDOMAINS)
  {
    setup("test/executable_mesh_extract_subdomains.xml");
    go();

    CHECK(S->HasMesh(""));
    CHECK(S->HasMesh("domain"));

    int has_upstream = 0;
    int set_size = S->GetMesh("domain")->getSetSize(
      "upstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    if (set_size > 0) {
      has_upstream += 1;
      CHECK(S->HasMesh("watershed:upstream"));
      CHECK_EQUAL(
        set_size,
        S->GetMesh("watershed:upstream")
          ->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    }

    int total_has_upstream = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &has_upstream, &total_has_upstream);
    CHECK(total_has_upstream > 0);

    // create a vector and fill
    std::map<std::string, Teuchos::RCP<MultiVector_type>> subdomain_vecs;
    auto ds = S->GetDomainSet("watershed");
    for (const auto& subdomain : *ds) {
      subdomain_vecs[subdomain] =
        Teuchos::rcp(new MultiVector_type(S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL,false), 1));
      if (subdomain == "watershed:upstream") {
        subdomain_vecs[subdomain]->putScalar(1.);
      } else {
        subdomain_vecs[subdomain]->putScalar(2.);
      }
    }

    // import to a global vector
    MultiVector_type vec(S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
    vec.putScalar(10.);

    for (const auto& subdomain : *ds) { ds->doImport(subdomain, *subdomain_vecs[subdomain], vec); }

    VerboseObject vo("EXTRACT_SUBDOMAINS", "extreme");
    std::cout << "WS:" << std::endl;
    subdomain_vecs["watershed:upstream"]->describe(*vo.os(), Teuchos::VERB_EXTREME);
    std::cout << "Domain:" << std::endl;
    vec.describe(*vo.os(), Teuchos::VERB_EXTREME);

    double result = vec.getVector(0)->normInf();
    CHECK_EQUAL(2.0, result); // not 10!
  }


#if 0
  TEST_FIXTURE(Runner, EXTRACT_SUBDOMAINS_SURFACE)
  {
    setup("test/executable_mesh_extract_subdomains_surface.xml");
    go();

    CHECK(S->HasMesh(""));
    CHECK(S->HasMesh("domain"));

    // check we got a valid surface mesh
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
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &has_surface, &total_has_surface);
    CHECK(total_has_surface > 0);

    // check we got upstream/downstream subdomains
    int has_upstream = 0;
    int set_size_us = S->GetMesh("domain")->getSetSize(
      "upstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    if (set_size_us > 0) {
      has_upstream += 1;
      CHECK(S->HasMesh("watershed:upstream"));
      CHECK_EQUAL(
        set_size_us,
        S->GetMesh("watershed:upstream")
          ->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    }

    int total_has_upstream = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &has_upstream, &total_has_upstream);
    CHECK(total_has_upstream > 0);

    // check we got upstream surface subdomain
    int has_surf_upstream = has_surface && has_upstream ? 1 : 0;
    if (has_surf_upstream) { CHECK(S->HasMesh("surface_watershed:upstream")); }
    int total_has_surf_upstream = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &has_surf_upstream, &total_has_surf_upstream);
    CHECK(total_has_surf_upstream > 0);

    // create a vector and fill
    std::map<std::string, Teuchos::RCP<MultiVector_type>> subdomain_vecs;
    auto ds = S->GetDomainSet("surface_watershed");
    for (const auto& subdomain : *ds) {
      subdomain_vecs[subdomain] =
        Teuchos::rcp(new MultiVector_type(S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL,false), 1));
      if (subdomain == "surface_watershed:upstream") {
        subdomain_vecs[subdomain]->putScalar(1.);
      } else {
        subdomain_vecs[subdomain]->putScalar(2.);
      }
    }

    // import to a global vector
    if (S->HasMesh("surface")) {
      MultiVector_type vec(S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
      vec.putScalar(10);

      for (const auto& subdomain : *ds) {
        ds->doImport(subdomain, *subdomain_vecs[subdomain], vec);
      }

      double result = vec.getVector(0)->normInf();
      CHECK_EQUAL(2.0, result); // not 10!

      result = vec.getVector(0)->norm2();
      CHECK_CLOSE(1.5, result, 1.e-6);
    }
  }


  TEST_FIXTURE(Runner, CONSTRUCT_COLUMNS)
  {
    setup("test/executable_mesh_construct_columns.xml");
    go();

    CHECK(S->HasMesh(""));
    CHECK(S->HasMesh("domain"));

    int has_surface = 0;
    int set_size = S->GetMesh("domain")->getSetSize(
      "surface", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    CHECK(set_size > 0); // columnar partitioned
    if (set_size > 0) {
      has_surface += 1;
      CHECK(S->HasMesh("surface"));
      CHECK_EQUAL(set_size,
                  S->GetMesh("surface")->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                      AmanziMesh::Parallel_kind::OWNED));
    }

    int num_cells = S->GetMesh("domain")->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                       AmanziMesh::Parallel_kind::OWNED);
    CHECK_EQUAL(0, num_cells % set_size);
    int ncells_per_column = num_cells / set_size;
    for (int col = 0; col != set_size; ++col) {
      // check that columns were made correctly
      CHECK_EQUAL(ncells_per_column, S->GetMesh("domain")->columns.getCells(col).size());
      CHECK_EQUAL(ncells_per_column + 1, S->GetMesh("domain")->columns.getFaces(col).size());

      // column mesh
      std::string col_name =
        Keys::getDomainInSet("column", S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL,false)->getGlobalElement(col));
      CHECK(S->HasMesh(col_name));
      CHECK_EQUAL(ncells_per_column,
                  S->GetMesh(col_name)->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                     AmanziMesh::Parallel_kind::OWNED));
    }

    {
      // construct col id vector
      MultiVector_type vec1(S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
      MultiVector_type vec2(S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL,false), 1);

      auto ds = S->GetDomainSet("column");
      int col = 0;
      for (const auto& subdomain : *ds) {
        // fill via import
        MultiVector_type vec_l(S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
        int index = Keys::getDomainSetIndex<int>(subdomain);
        vec_l.putScalar((double)index);
        ds->doImport(subdomain, vec_l, vec2);

        // fill via column
        for (const auto& c : S->GetMesh("domain")->columns.getCells(col)) {
          vec1.replaceLocalValue(c, 0, index);
        }
        col++;
      }

      // check they are the same
      vec1.update(-1, vec2, 1);
      double norm = vec1.getVector(0)->normInf();
      CHECK_CLOSE(0., norm, 1.e-10);
    }

    // check that column surfaces were made correctly
    for (int col = 0; col != set_size; ++col) {
      // column mesh
      std::string surf_col_name =
        Keys::getDomainInSet("surface_column", S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL,false)->getGlobalElement(col));
      CHECK(S->HasMesh(surf_col_name));
      CHECK_EQUAL(
        1,
        S->GetMesh(surf_col_name)
          ->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    }


    {
      // construct col id vector
      MultiVector_type vec1(S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
      MultiVector_type vec2(S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL,false), 1);

      auto ds = S->GetDomainSet("surface_column");
      int col = 0;
      for (const auto& subdomain : *ds) {
        // fill via import
        MultiVector_type vec_l(S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL,false), 1);
        int index = Keys::getDomainSetIndex<int>(subdomain);
        vec_l.putScalar((double)index);
        ds->doImport(subdomain, vec_l, vec2);

        // fill via column
        vec1.replaceLocalValue(col, 0, index);
        col++;
      }

      // check they are the same
      vec1.update(-1, vec2, 1);
      double norm = vec1.getVector(0)->normInf();
      CHECK_CLOSE(0., norm, 1.e-10);
    }
  }
#endif
}
