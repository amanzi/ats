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
#include "Reductions.hh"
#include "ats_mesh_factory.hh"

using namespace Amanzi;

void
parallel_print(const Comm_ptr_type& comm, const std::string& str)
{
  for (int i = 0; i != comm->getSize(); ++i) {
    if (i == comm->getRank()) {
      std::cout << "[" << i << "/" << comm->getSize() << "] " << str;
      std::flush(std::cout);
    }
    comm->barrier();
  }
}


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
    S = Teuchos::rcp(new State(Teuchos::sublist(plist, "state")));
  }

  void go() { ATS::Mesh::createMeshes(plist, comm, gm, *S); }

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

    // validate the mesh?
    auto& domain = *S->GetMesh("domain");
    auto ncells = domain.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    std::stringstream str;
    str << "Created a DOMAIN mesh with " << ncells << " cells" << std::endl;
    parallel_print(comm, str.str());

    int ncells_g;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ncells, &ncells_g);
    CHECK_EQUAL(8400, ncells_g);

    // validate the upstream set
    auto upstream_count = domain.getSetSize(
      "upstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    std::stringstream str2;
    str2 << "Upstream Region has " << upstream_count << " cells" << std::endl;
    parallel_print(comm, str2.str());
    int upstream_count_g;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &upstream_count, &upstream_count_g);
    CHECK_EQUAL(4200, upstream_count_g);

    // validate the downstream set
    auto downstream_count = domain.getSetSize(
      "downstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    std::stringstream str3;
    str3 << "Downstream Region has " << downstream_count << " cells" << std::endl;
    parallel_print(comm, str3.str());
    int downstream_count_g;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &downstream_count, &downstream_count_g);
    CHECK_EQUAL(4200, downstream_count_g);

    // validate hte upstream mesh
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
    } else {
      CHECK(!S->HasMesh("watershed:upstream"));
    }

    int total_has_upstream = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &has_upstream, &total_has_upstream);
    CHECK(total_has_upstream > 0);

    // NOTE: this is NOT gauranteed for all partitionings, so this may fail
    // eventually.  This check simply makes sure the test is valid and hlepful
    // by having at least one rank without the upstream domain
    if (comm->getSize() > 1) CHECK(total_has_upstream < comm->getSize());
    if (comm->getRank() == 0)
      std::cout << "Upstream found on " << total_has_upstream << " / " << comm->getSize()
                << " ranks" << std::endl;

    // check the right sizes of the upstream mesh
    if (has_upstream) {
      auto& up_mesh = *S->GetMesh("watershed:upstream");
      int ncells_upstream = up_mesh.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
      int ncells_upstream_g;
      Teuchos::reduceAll(
        *up_mesh.getComm(), Teuchos::REDUCE_SUM, 1, &ncells_upstream, &ncells_upstream_g);
      CHECK_EQUAL(4200, ncells_upstream_g);
    }

    // create a vector and fill by domain set
    std::map<std::string, Teuchos::RCP<MultiVector_type>> subdomain_vecs;
    auto ds = S->GetDomainSet("watershed");
    for (const auto& subdomain : *ds) {
      subdomain_vecs[subdomain] = Teuchos::rcp(new MultiVector_type(
        S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL, false), 1));
      if (subdomain == "watershed:upstream") {
        subdomain_vecs[subdomain]->putScalar(1.);
      } else {
        subdomain_vecs[subdomain]->putScalar(2.);
      }
    }

    // import to a global vector
    MultiVector_type vec(S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL, false), 1);
    vec.putScalar(-1);

    for (const auto& subdomain : *ds) { ds->doImport(subdomain, *subdomain_vecs[subdomain], vec); }

    auto result = Reductions::reduceAllMinLoc(*vec.getVector(0));
    CHECK_EQUAL(1.0, result.val);

    double result2 = vec.getVector(0)->norm1();
    CHECK_CLOSE(4200 * 1 + 4200 * 2, result2, 1.e-6);
  }

#if 0
  TEST_FIXTURE(Runner, EXTRACT_SUBDOMAINS_SURFACE)
  {
    setup("test/executable_mesh_extract_subdomains_surface.xml");
    go();

    CHECK(S->HasMesh(""));
    CHECK(S->HasMesh("domain"));

    // validate the mesh?
    auto& domain = *S->GetMesh("domain");
    auto ncells = domain.ncells_owned;
    std::stringstream str;
    str << "Created a DOMAIN mesh with " << ncells << " cells" << std::endl;
    parallel_print(comm, str.str());

    int ncells_g;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &ncells, &ncells_g);
    CHECK_EQUAL(8400, ncells_g);

    // validate the surface mesh?
    int has_surface = 0;
    int set_size = S->GetMesh("domain")->getSetSize(
      "surface", AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
    if (set_size > 0) {
      has_surface += 1;
      CHECK(S->HasMesh("surface"));
      CHECK_EQUAL(set_size,
                  S->GetMesh("surface")->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                        AmanziMesh::Parallel_kind::OWNED));
    } else {
      CHECK(!S->HasMesh("surface"));
    }
    int total_has_surface = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &has_surface, &total_has_surface);
    CHECK(total_has_surface > 0);
    // horizontal decomposition means all ranks have surface
    CHECK(total_has_surface == comm->getSize());

    // validate the upstream set
    auto upstream_count = domain.getSetSize(
      "upstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    std::stringstream str2;
    str2 << "Upstream Region has " << upstream_count << " cells" << std::endl;
    parallel_print(comm, str2.str());
    int upstream_count_g;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &upstream_count, &upstream_count_g);
    CHECK_EQUAL(4200, upstream_count_g);

    // validate the downstream set
    auto downstream_count = domain.getSetSize(
      "downstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    std::stringstream str3;
    str3 << "Downstream Region has " << downstream_count << " cells" << std::endl;
    parallel_print(comm, str3.str());
    int downstream_count_g;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &downstream_count, &downstream_count_g);
    CHECK_EQUAL(4200, downstream_count_g);

    // validate the upstream mesh
    int has_upstream = 0;
    int us_set_size = S->GetMesh("domain")->getSetSize(
      "upstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    if (us_set_size > 0) {
      has_upstream += 1;
      CHECK(S->HasMesh("watershed:upstream"));
      CHECK_EQUAL(
        us_set_size,
        S->GetMesh("watershed:upstream")
          ->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    } else {
      CHECK(!S->HasMesh("watershed:upstream"));
    }

    int total_has_upstream = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &has_upstream, &total_has_upstream);
    CHECK(total_has_upstream > 0);

    // NOTE: this is NOT gauranteed for all partitionings, so this may fail
    // eventually.  This check simply makes sure the test is valid and hlepful
    // by having at least one rank without the upstream domain
    if (comm->getSize() > 1) CHECK(total_has_upstream < comm->getSize());
    if (comm->getRank() == 0)
      std::cout << "Upstream found on " << total_has_upstream << " / " << comm->getSize()
                << " ranks" << std::endl;

    // validate the upstream
    int has_surf_upstream = has_surface && has_upstream ? 1 : 0;
    if (has_surf_upstream) {
      CHECK(S->HasMesh("surface_watershed:upstream"));
    } else {
      CHECK(!S->HasMesh("surface_watershed:upstream"));
    }
    int total_has_surf_upstream = 0;
    Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, 1, &has_surf_upstream, &total_has_surf_upstream);
    CHECK(total_has_surf_upstream > 0);
    CHECK_EQUAL(total_has_upstream, total_has_surf_upstream); // horizontal decomposition

    // create a vector and fill
    std::map<std::string, Teuchos::RCP<MultiVector_type>> subdomain_vecs;
    auto ds = S->GetDomainSet("surface_watershed");
    for (const auto& subdomain : *ds) {
      subdomain_vecs[subdomain] = Teuchos::rcp(new MultiVector_type(
        S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL, false), 1));
      if (subdomain == "surface_watershed:upstream") {
        subdomain_vecs[subdomain]->putScalar(1.);
      } else {
        subdomain_vecs[subdomain]->putScalar(2.);
      }
    }

    // import to a global vector
    if (S->HasMesh("surface")) {
      MultiVector_type vec(S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false),
                             1);
      vec.putScalar(-1);

      for (const auto& subdomain : *ds) {
        ds->doImport(subdomain, *subdomain_vecs[subdomain], vec);
      }

      auto result = Reductions::reduceAllMinLoc(*vec.getVector(0));
      CHECK_EQUAL(1.0, result.val);

      double result2 = vec.getVector(0)->norm1();
      CHECK_CLOSE(1 * 420 + 2 * 420, result2, 1.e-6);
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
      CHECK_EQUAL(ncells_per_column, S->GetMesh("domain")->columns->getCells<MemSpace_kind::HOST>(col).size());
      CHECK_EQUAL(ncells_per_column + 1, S->GetMesh("domain")->columns->getFaces<MemSpace_kind::HOST>(col).size());

      // column mesh
      std::string col_name = Keys::getDomainInSet(
        "column", S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false)->getGlobalElement(col));
      CHECK(S->HasMesh(col_name));
      CHECK_EQUAL(ncells_per_column,
                  S->GetMesh(col_name)->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                       AmanziMesh::Parallel_kind::OWNED));
    }

    {
      // construct col id vector
      MultiVector_type vec1(S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL, false),
                              1);
      MultiVector_type vec2(S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL, false),
                              1);

      auto ds = S->GetDomainSet("column");
      int col = 0;
      for (const auto& subdomain : *ds) {
        // fill via import
        MultiVector_type vec_l(
          S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL, false), 1);
        int index = Keys::getDomainSetIndex<int>(subdomain);
        vec_l.putScalar((double)index);
        ds->doImport(subdomain, vec_l, vec2);

        // fill via column
        for (const auto& c : S->GetMesh("domain")->columns->getCells<MemSpace_kind::HOST>(col)) {
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
      std::string surf_col_name = Keys::getDomainInSet(
        "surface_column",
        S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false)->getGlobalElement(col));
      CHECK(S->HasMesh(surf_col_name));
      CHECK_EQUAL(
        1,
        S->GetMesh(surf_col_name)
          ->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    }


    {
      // construct col id vector
      MultiVector_type vec1(S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false),
                              1);
      MultiVector_type vec2(S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false),
                              1);

      auto ds = S->GetDomainSet("surface_column");
      int col = 0;
      for (const auto& subdomain : *ds) {
        // fill via import
        MultiVector_type vec_l(
          S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL, false), 1);
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
