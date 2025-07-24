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

void
parallel_print(const Comm_ptr_type& comm, const std::string& str)
{
  for (int i = 0; i != comm->NumProc(); ++i) {
    if (i == comm->MyPID()) {
      std::cout << "[" << i << "/" << comm->NumProc() << "] " << str;
      std::flush(std::cout);
    }
    comm->Barrier();
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
    S = Teuchos::rcp(new State(plist->sublist("state")));
  }

  void go() { ATS::Mesh::createMeshes(plist, comm, gm, *S); }

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

  TEST_FIXTURE(Runner, EXTRACT_SUBDOMAINS)
  {
    setup("test/executable_mesh_extract_subdomains.xml");
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
    comm->SumAll(&ncells, &ncells_g, 1);
    CHECK_EQUAL(8400, ncells_g);

    // validate the upstream set
    auto upstream_count = domain.getSetSize(
      "upstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    std::stringstream str2;
    str2 << "Upstream Region has " << upstream_count << " cells" << std::endl;
    parallel_print(comm, str2.str());
    int upstream_count_g;
    comm->SumAll(&upstream_count, &upstream_count_g, 1);
    CHECK_EQUAL(4200, upstream_count_g);

    // validate the downstream set
    auto downstream_count = domain.getSetSize(
      "downstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    std::stringstream str3;
    str3 << "Downstream Region has " << downstream_count << " cells" << std::endl;
    parallel_print(comm, str3.str());
    int downstream_count_g;
    comm->SumAll(&downstream_count, &downstream_count_g, 1);
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
    comm->SumAll(&has_upstream, &total_has_upstream, 1);
    CHECK(total_has_upstream > 0);

    // NOTE: this is NOT gauranteed for all partitionings, so this may fail
    // eventually.  This check simply makes sure the test is valid and hlepful
    // by having at least one rank without the upstream domain
    if (comm->NumProc() > 1) CHECK(total_has_upstream < comm->NumProc());
    if (comm->MyPID() == 0)
      std::cout << "Upstream found on " << total_has_upstream << " / " << comm->NumProc()
                << " ranks" << std::endl;

    // check the right sizes of the upstream mesh
    if (has_upstream) {
      auto& up_mesh = *S->GetMesh("watershed:upstream");
      int ncells_upstream = up_mesh.ncells_owned;
      int ncells_upstream_g;
      up_mesh.getComm()->SumAll(&ncells_upstream, &ncells_upstream_g, 1);
      CHECK_EQUAL(4200, ncells_upstream_g);
    }

    // create a vector and fill by domain set
    std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> subdomain_vecs;
    auto ds = S->GetDomainSet("watershed");
    for (const auto& subdomain : *ds) {
      subdomain_vecs[subdomain] = Teuchos::rcp(new Epetra_MultiVector(
        S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL, false), 1));
      if (subdomain == "watershed:upstream") {
        subdomain_vecs[subdomain]->PutScalar(1.);
      } else {
        subdomain_vecs[subdomain]->PutScalar(2.);
      }
    }

    // import to a global vector
    Epetra_MultiVector vec(S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL, false), 1);
    vec.PutScalar(-1);

    for (const auto& subdomain : *ds) {
      ds->doImport(subdomain, *subdomain_vecs[subdomain], vec);
    }

    double result;
    vec.MinValue(&result);
    CHECK_EQUAL(1.0, result);

    vec.Norm1(&result);
    CHECK_CLOSE(4200 * 1 + 4200 * 2, result, 1.e-6);
  }


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
    comm->SumAll(&ncells, &ncells_g, 1);
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
    comm->SumAll(&has_surface, &total_has_surface, 1);
    CHECK(total_has_surface > 0);
    // horizontal decomposition means all ranks have surface
    CHECK(total_has_surface == comm->NumProc());

    // validate the upstream set
    auto upstream_count = domain.getSetSize(
      "upstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    std::stringstream str2;
    str2 << "Upstream Region has " << upstream_count << " cells" << std::endl;
    parallel_print(comm, str2.str());
    int upstream_count_g;
    comm->SumAll(&upstream_count, &upstream_count_g, 1);
    CHECK_EQUAL(4200, upstream_count_g);

    // validate the downstream set
    auto downstream_count = domain.getSetSize(
      "downstream", AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    std::stringstream str3;
    str3 << "Downstream Region has " << downstream_count << " cells" << std::endl;
    parallel_print(comm, str3.str());
    int downstream_count_g;
    comm->SumAll(&downstream_count, &downstream_count_g, 1);
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
    comm->SumAll(&has_upstream, &total_has_upstream, 1);
    CHECK(total_has_upstream > 0);

    // NOTE: this is NOT gauranteed for all partitionings, so this may fail
    // eventually.  This check simply makes sure the test is valid and hlepful
    // by having at least one rank without the upstream domain
    if (comm->NumProc() > 1) CHECK(total_has_upstream < comm->NumProc());
    if (comm->MyPID() == 0)
      std::cout << "Upstream found on " << total_has_upstream << " / " << comm->NumProc()
                << " ranks" << std::endl;

    // validate the upstream
    int has_surf_upstream = has_surface && has_upstream ? 1 : 0;
    if (has_surf_upstream) {
      CHECK(S->HasMesh("surface_watershed:upstream"));
    } else {
      CHECK(!S->HasMesh("surface_watershed:upstream"));
    }
    int total_has_surf_upstream = 0;
    comm->SumAll(&has_surf_upstream, &total_has_surf_upstream, 1);
    CHECK(total_has_surf_upstream > 0);
    CHECK_EQUAL(total_has_upstream, total_has_surf_upstream); // horizontal decomposition

    // create a vector and fill
    std::map<std::string, Teuchos::RCP<Epetra_MultiVector>> subdomain_vecs;
    auto ds = S->GetDomainSet("surface_watershed");
    for (const auto& subdomain : *ds) {
      subdomain_vecs[subdomain] = Teuchos::rcp(new Epetra_MultiVector(
        S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL, false), 1));
      if (subdomain == "surface_watershed:upstream") {
        subdomain_vecs[subdomain]->PutScalar(1.);
      } else {
        subdomain_vecs[subdomain]->PutScalar(2.);
      }
    }

    // import to a global vector
    if (S->HasMesh("surface")) {
      Epetra_MultiVector vec(S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false),
                             1);
      vec.PutScalar(-1);

      for (const auto& subdomain : *ds) {
        ds->doImport(subdomain, *subdomain_vecs[subdomain], vec);
      }

      double result;
      vec.MinValue(&result);
      CHECK_EQUAL(1.0, result);

      vec.Norm1(&result);
      CHECK_CLOSE(1 * 420 + 2 * 420, result, 1.e-6);
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
      std::string col_name = Keys::getDomainInSet(
        "column", S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false).GID(col));
      CHECK(S->HasMesh(col_name));
      CHECK_EQUAL(ncells_per_column,
                  S->GetMesh(col_name)->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                       AmanziMesh::Parallel_kind::OWNED));
    }

    {
      // construct col id vector
      Epetra_MultiVector vec1(S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL, false),
                              1);
      Epetra_MultiVector vec2(S->GetMesh("domain")->getMap(AmanziMesh::Entity_kind::CELL, false),
                              1);

      auto ds = S->GetDomainSet("column");
      int col = 0;
      for (const auto& subdomain : *ds) {
        // fill via import
        Epetra_MultiVector vec_l(
          S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL, false), 1);
        int index = Keys::getDomainSetIndex<int>(subdomain);
        vec_l.PutScalar((double)index);
        ds->doImport(subdomain, vec_l, vec2);

        // fill via column
        for (const auto& c : S->GetMesh("domain")->columns.getCells(col)) {
          vec1[0][c] = index;
        }
        col++;
      }

      // check they are the same
      vec1.Update(-1, vec2, 1);
      double norm;
      vec1.NormInf(&norm);
      CHECK_CLOSE(0., norm, 1.e-10);
    }

    // check that column surfaces were made correctly
    for (int col = 0; col != set_size; ++col) {
      // column mesh
      std::string surf_col_name = Keys::getDomainInSet(
        "surface_column",
        S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false).GID(col));
      CHECK(S->HasMesh(surf_col_name));
      CHECK_EQUAL(
        1,
        S->GetMesh(surf_col_name)
          ->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED));
    }


    {
      // construct col id vector
      Epetra_MultiVector vec1(S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false),
                              1);
      Epetra_MultiVector vec2(S->GetMesh("surface")->getMap(AmanziMesh::Entity_kind::CELL, false),
                              1);

      auto ds = S->GetDomainSet("surface_column");
      int col = 0;
      for (const auto& subdomain : *ds) {
        // fill via import
        Epetra_MultiVector vec_l(
          S->GetMesh(subdomain)->getMap(AmanziMesh::Entity_kind::CELL, false), 1);
        int index = Keys::getDomainSetIndex<int>(subdomain);
        vec_l.PutScalar((double)index);
        ds->doImport(subdomain, vec_l, vec2);

        // fill via column
        vec1[0][col] = index;
        col++;
      }

      // check they are the same
      vec1.Update(-1, vec2, 1);
      double norm;
      vec1.NormInf(&norm);
      CHECK_CLOSE(0., norm, 1.e-10);
    }
  }
}
