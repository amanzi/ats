/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*

*/

#include <iostream>
#include <string>
#include <vector>

#include "AmanziComm.hh"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "UnitTest++.h"

// Amanzi
#include "exceptions.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "EvaluatorSecondary.hh"

#include "mpc_surface_subsurface_helpers.hh"
#include "mpc_coupled_water.hh"
#include "ats_mesh_factory.hh"
#include "richards.hh"
#include "overland_pressure.hh"

using namespace Amanzi;
using namespace ATS;
using namespace Amanzi::AmanziMesh;
using namespace Amanzi::AmanziGeometry;
using namespace Amanzi::ATS_Physics;

struct CoupledWaterProblem {
 public:
  Teuchos::RCP<Teuchos::ParameterList> plist;
  Teuchos::RCP<State> S;
  Teuchos::RCP<MPCCoupledWater> pk;
  Teuchos::RCP<Flow::Richards> pk_richards;
  Teuchos::RCP<Flow::OverlandPressureFlow> pk_overland;
  Teuchos::RCP<TreeVector> soln, soln2, soln_old, res, res2, dsoln;
  Comm_ptr_type comm;

  CoupledWaterProblem() {}

  void init()
  {
    comm = getDefaultComm();

    // create state.
    S = Teuchos::rcp(new State(plist->sublist("state")));
    soln = Teuchos::rcp(new TreeVector(comm));

    // create meshes
    auto& regions_list = plist->sublist("regions");
    auto gm = Teuchos::rcp(new GeometricModel(3, regions_list, *comm));
    ATS::Mesh::createMeshes(plist, comm, gm, *S);

    // create the PK
    Teuchos::ParameterList pk_tree_list("PK tree");
    pk_tree_list.sublist("coupled water").set("PK type", "coupled water");
    pk_tree_list.sublist("coupled water")
      .sublist("subsurface flow")
      .set("PK type", "richards flow");
    pk_tree_list.sublist("coupled water")
      .sublist("surface flow")
      .set("PK type", "overland flow, pressure basis");

    Amanzi::PKFactory pk_factory;
    auto pk_as_pk = pk_factory.CreatePK("coupled water", pk_tree_list, plist, S, soln);
    pk = Teuchos::rcp_dynamic_cast<MPCCoupledWater>(pk_as_pk);
    pk->parseParameterList();
    AMANZI_ASSERT(pk.get());

    // setup stage
    // common constants
    S->Require<double>("atmospheric_pressure", Tags::DEFAULT, "coordinator");
    S->Require<AmanziGeometry::Point>("gravity", Tags::DEFAULT, "coordinator");

    // needed other times
    S->require_time(Tags::CURRENT);
    S->require_time(Tags::NEXT);

    pk->set_tags(Amanzi::Tags::CURRENT, Amanzi::Tags::NEXT);
    pk->Setup();
    S->Setup();

    // initialization stage
    S->set_time(Amanzi::Tags::CURRENT, 0.);
    S->set_time(Amanzi::Tags::NEXT, 0.);
    S->set_cycle(0);
    S->InitializeFields();
    pk->Initialize();
    pk->CommitStep(0., 0., Tags::NEXT);

    S->InitializeEvaluators();
    S->InitializeFieldCopies();
    S->CheckAllFieldsInitialized();
    pk->CommitStep(0, 0, Tags::NEXT);
    if (S->get_cycle() == -1) S->advance_cycle();

    S->set_time(Amanzi::Tags::NEXT, 8640.);

    // set up other needed vectors
    pk->State_to_Solution(Amanzi::Tags::NEXT, *soln);
    soln_old = Teuchos::rcp(new TreeVector(*soln));
    pk->State_to_Solution(Amanzi::Tags::CURRENT, *soln_old);

    soln2 = Teuchos::rcp(new TreeVector(*soln));
    *soln2 = *soln;

    dsoln = Teuchos::rcp(new TreeVector(*soln));
    dsoln->PutScalar(0.);

    res = Teuchos::rcp(new TreeVector(*soln));
    res->PutScalar(0.);

    res2 = Teuchos::rcp(new TreeVector(*soln));
    res2->PutScalar(0.);

    // stash pointers to the sub-pks
    pk_richards = Teuchos::rcp_dynamic_cast<Flow::Richards>(pk->get_subpk(0));
    pk_overland = Teuchos::rcp_dynamic_cast<Flow::OverlandPressureFlow>(pk->get_subpk(1));
  }
};


SUITE(EXECUTABLE_COUPLED_WATER)
{
  TEST_FIXTURE(CoupledWaterProblem, CONSTRUCT)
  {
    std::cout << std::endl
              << std::endl
              << "CoupledWaterProblem Construction" << std::endl
              << "============================================================" << std::endl;
    plist = Teuchos::getParametersFromXmlFile("test/executable_coupled_water1.xml");
    init();
  };


  //
  // This test fixes the viscosity of water at very large, effectively turning
  // off Darcy flux and thereby testing the preconditioner's accumulation term
  // only.
  //
  TEST_FIXTURE(CoupledWaterProblem, PRECONDITIONER1_ACCUMULATION)
  {
    std::cout << std::endl
              << std::endl
              << "CoupledWaterProblem Precon1: Accumulation" << std::endl
              << "============================================================" << std::endl;
    double eps = .001;
    Entity_ID c = 1;

    plist = Teuchos::getParametersFromXmlFile("test/executable_coupled_water1.xml");
    Teuchos::Array<int> dc(1);
    dc[0] = c;
    plist->sublist("PKs").sublist("subsurface flow").set<Teuchos::Array<int>>("debug cells", dc);
    init();


    // evaluate residual
    std::cout << "EVALUATING RES: BASE" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res);
    pk->ErrorNorm(soln, res);

    // update the preconditioner
    pk->UpdatePreconditioner(8640, soln, 8640);

    // Now perturb and eval residual again
    // -- PC1: change a cell value
    CHECK_EQUAL(5 * 25, dsoln->SubVector(0)->Data()->ViewComponent("cell", false)->MyLength());
    (*dsoln->SubVector(0)->Data()->ViewComponent("cell", false))[0][c] = eps;
    soln->Update(1, *dsoln, 1);

    // -- mark pressure as changed
    pk->ChangedSolutionPK(Tags::NEXT);

    // -- evaluate the residual again
    std::cout << std::endl
              << "EVALUATING RES: PERTURBED" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res2);
    pk->ErrorNorm(soln, res2);

    // Compute r(p + dp) - r(p) ~= J * dp
    // -- diff residuals
    res2->Update(-1, *res, 1.);

    // -- apply the PC
    std::cout << std::endl
              << "APPLY PC" << std::endl
              << "--------------------------------------------------" << std::endl;
    res->PutScalar(0.);
    pk->preconditioner()->Apply(*dsoln->SubVector(0)->Data(), *res->SubVector(0)->Data());
    CopySubsurfaceToSurface(*res->SubVector(0)->Data(), *res->SubVector(1)->Data());


    pk_richards->debugger()->WriteVector("dres", res2->SubVector(0)->Data().ptr());
    pk_richards->debugger()->WriteVector("Jdp", res->SubVector(0)->Data().ptr());

    // -- subtract J dp
    double norm;
    res->Update(-1 / eps, *res2, 1. / eps);
    res->NormInf(&norm);
    CHECK_CLOSE(0., norm, 1.e-3);
  };


  //
  // This test fixes the water content at 1, thereby turning off accumulation and
  // testing diffusion only, with no Jacobian terms and no gravity.
  //
  TEST_FIXTURE(CoupledWaterProblem, PRECONDITIONER2_DIFFUSION_NO_GRAVITY)
  {
    std::cout << std::endl
              << std::endl
              << "CoupledWaterProblem Precon2: Diffusion 0 Gravity" << std::endl
              << "============================================================" << std::endl;
    double eps = .0001;
    Entity_ID c = 1;

    plist = Teuchos::getParametersFromXmlFile("test/executable_coupled_water2.xml");
    Teuchos::Array<int> dc(1);
    dc[0] = c;
    plist->sublist("PKs").sublist("subsurface flow").set<Teuchos::Array<int>>("debug cells", dc);
    init();

    // evaluate residual
    std::cout << "EVALUATING RES: BASE" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res);
    pk->ErrorNorm(soln, res);

    // update the preconditioner
    pk->UpdatePreconditioner(8640, soln, 8640);

    // Now perturb and eval residual again
    // -- PC1: change a cell value
    CHECK_EQUAL(5 * 25, dsoln->SubVector(0)->Data()->ViewComponent("cell", false)->MyLength());
    (*dsoln->SubVector(0)->Data()->ViewComponent("cell", false))[0][c] = eps;
    soln->Update(1, *dsoln, 1);

    // -- mark pressure as changed
    pk->ChangedSolutionPK(Tags::NEXT);

    // -- fix krel to not include the effect of dkr/dp
    pk_richards->set_fixed_kr();

    // -- evaluate the residual again
    std::cout << std::endl
              << "EVALUATING RES: PERTURBED" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res2);
    pk->ErrorNorm(soln, res2);

    // Compute r(p + dp) - r(p) ~= J * dp
    // -- diff residuals
    res2->Update(-1, *res, 1.);

    // -- apply the PC
    std::cout << std::endl
              << "APPLY PC" << std::endl
              << "--------------------------------------------------" << std::endl;
    res->PutScalar(0.);
    pk->preconditioner()->Apply(*dsoln->SubVector(0)->Data(), *res->SubVector(0)->Data());
    CopySubsurfaceToSurface(*res->SubVector(0)->Data(), *res->SubVector(1)->Data());


    pk_richards->debugger()->WriteVector("dres", res2->SubVector(0)->Data().ptr(), true);
    pk_richards->debugger()->WriteVector("Jdp", res->SubVector(0)->Data().ptr(), true);

    // -- subtract J dp
    double norm;
    res->Update(-1 / eps, *res2, 1. / eps);
    res->NormInf(&norm);
    CHECK_CLOSE(0., norm, 1.e-3);
  };


  //
  // This test fixes the water content at 1, thereby turning off accumulation and
  // testing diffusion only, with no Jacobian terms and no gravity.
  //
  TEST_FIXTURE(CoupledWaterProblem, PRECONDITIONER3_DIFFUSION_WITH_GRAVITY)
  {
    std::cout << std::endl
              << std::endl
              << "CoupledWaterProblem Precon3: Diffusion With Gravity" << std::endl
              << "============================================================" << std::endl;
    double eps = .0001;
    Entity_ID c = 1;

    plist = Teuchos::getParametersFromXmlFile("test/executable_coupled_water3.xml");
    Teuchos::Array<int> dc(1);
    dc[0] = c;
    plist->sublist("PKs").sublist("subsurface flow").set<Teuchos::Array<int>>("debug cells", dc);
    init();

    auto db = pk_richards->debugger();

    // evaluate residual
    std::cout << "EVALUATING RES: BASE" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res);
    db->WriteVector("mat_diff rhs", pk_richards->get_operator()->rhs().ptr(), true);
    pk->ErrorNorm(soln, res);


    // update the preconditioner
    pk->UpdatePreconditioner(8640, soln, 8640);

    // Now perturb and eval residual again
    // -- PC1: change a cell value
    CHECK_EQUAL(5 * 25, dsoln->SubVector(0)->Data()->ViewComponent("cell", false)->MyLength());
    (*dsoln->SubVector(0)->Data()->ViewComponent("cell", false))[0][c] = eps;
    soln->Update(1, *dsoln, 1);

    // -- mark pressure as changed
    pk->ChangedSolutionPK(Tags::NEXT);

    // -- fix krel to not include the effect of dkr/dp
    pk_richards->set_fixed_kr();

    // -- evaluate the residual again
    std::cout << std::endl
              << "EVALUATING RES: PERTURBED" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res2);
    db->WriteVector("mat_diff rhs", pk_richards->get_operator()->rhs().ptr(), true);
    pk->ErrorNorm(soln, res2);

    // Compute r(p + dp) - r(p) ~= J * dp
    // -- diff residuals
    res2->Update(-1, *res, 1.);

    // -- apply the PC
    std::cout << std::endl
              << "APPLY PC" << std::endl
              << "--------------------------------------------------" << std::endl;
    res->PutScalar(0.);
    pk->preconditioner()->Apply(*dsoln->SubVector(0)->Data(), *res->SubVector(0)->Data());
    CopySubsurfaceToSurface(*res->SubVector(0)->Data(), *res->SubVector(1)->Data());


    db->WriteVector("dres", res2->SubVector(0)->Data().ptr(), true);
    db->WriteVector("Jdp", res->SubVector(0)->Data().ptr(), true);

    // -- subtract J dp
    double norm;
    res->Update(-1 / eps, *res2, 1. / eps);
    res->NormInf(&norm);
    CHECK_CLOSE(0., norm, 1.e-3);
  };


  //
  // This test fixes the water content at 1, thereby turning off accumulation and
  // testing diffusion only, with no Jacobian terms and no gravity.
  //
  TEST_FIXTURE(CoupledWaterProblem, PRECONDITIONER4_DIFFUSION_WITH_GRAVITY_SURF)
  {
    std::cout << std::endl
              << std::endl
              << "CoupledWaterProblem Precon3: Diffusion With Gravity (CELL at top of domain)"
              << std::endl
              << "============================================================" << std::endl;
    double eps = .0001;
    Entity_ID c = 0;

    plist = Teuchos::getParametersFromXmlFile("test/executable_coupled_water3.xml");
    Teuchos::Array<int> dc(1);
    dc[0] = c;
    plist->sublist("PKs").sublist("subsurface flow").set<Teuchos::Array<int>>("debug cells", dc);
    init();

    auto db = pk_richards->debugger();

    // evaluate residual
    std::cout << "EVALUATING RES: BASE" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res);
    db->WriteVector("mat_diff rhs", pk_richards->get_operator()->rhs().ptr(), true);
    pk->ErrorNorm(soln, res);


    // update the preconditioner
    pk->UpdatePreconditioner(8640, soln, 8640);

    // Now perturb and eval residual again
    // -- PC1: change a cell value
    CHECK_EQUAL(5 * 25, dsoln->SubVector(0)->Data()->ViewComponent("cell", false)->MyLength());
    (*dsoln->SubVector(0)->Data()->ViewComponent("cell", false))[0][c] = eps;
    soln->Update(1, *dsoln, 1);

    // -- mark pressure as changed
    pk->ChangedSolutionPK(Tags::NEXT);

    // -- fix krel to not include the effect of dkr/dp
    pk_richards->set_fixed_kr();

    // -- evaluate the residual again
    std::cout << std::endl
              << "EVALUATING RES: PERTURBED" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res2);
    db->WriteVector("mat_diff rhs", pk_richards->get_operator()->rhs().ptr(), true);
    pk->ErrorNorm(soln, res2);

    // Compute r(p + dp) - r(p) ~= J * dp
    // -- diff residuals
    res2->Update(-1, *res, 1.);

    // -- apply the PC
    std::cout << std::endl
              << "APPLY PC" << std::endl
              << "--------------------------------------------------" << std::endl;
    res->PutScalar(0.);
    pk->preconditioner()->Apply(*dsoln->SubVector(0)->Data(), *res->SubVector(0)->Data());

    db->WriteVector("dres", res2->SubVector(0)->Data().ptr(), true);
    db->WriteVector("Jdp", res->SubVector(0)->Data().ptr(), true);

    // -- subtract J dp
    double norm;
    res->Update(-1 / eps, *res2, 1. / eps);
    res->NormInf(&norm);
    CHECK_CLOSE(0., norm, 1.e-3);
    if (norm > 1.e-3) {
      std::cout << "res(p + dp) - res(p) - J(dp):" << std::endl;
      res->Print(std::cout);
    }
  };

  //
  // This test fixes the water content at 1, thereby turning off accumulation and
  // testing diffusion only, with no Jacobian terms and no gravity.
  //
  TEST_FIXTURE(CoupledWaterProblem, PRECONDITIONER4_DIFFUSION_WITH_GRAVITY_SURF_FACE)
  {
    std::cout << std::endl
              << std::endl
              << "CoupledWaterProblem Precon3: Diffusion With Gravity (surface face)" << std::endl
              << "============================================================" << std::endl;
    double eps = .0001;
    Entity_ID f = 0;

    plist = Teuchos::getParametersFromXmlFile("test/executable_coupled_water3.xml");
    Teuchos::Array<int> df(1);
    df[0] = f;
    plist->sublist("PKs").sublist("subsurface flow").set<Teuchos::Array<int>>("debug faces", df);
    init();

    auto db = pk_richards->debugger();

    // evaluate residual
    std::cout << "EVALUATING RES: BASE" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res);
    db->WriteVector("mat_diff rhs", pk_richards->get_operator()->rhs().ptr(), true);
    pk->ErrorNorm(soln, res);


    // update the preconditioner
    pk->UpdatePreconditioner(8640, soln, 8640);

    // Now perturb and eval residual again
    // -- PC1: change a cell value
    CHECK_EQUAL(5 * 25, dsoln->SubVector(0)->Data()->ViewComponent("cell", false)->MyLength());
    (*dsoln->SubVector(0)->Data()->ViewComponent("face", false))[0][f] = eps;
    soln->Update(1, *dsoln, 1);
    CopySubsurfaceToSurface(*soln->SubVector(0)->Data(), *soln->SubVector(1)->Data());

    // -- mark pressure as changed
    pk->ChangedSolutionPK(Tags::NEXT);

    // -- fix krel to not include the effect of dkr/dp
    pk_richards->set_fixed_kr();

    // -- evaluate the residual again
    std::cout << std::endl
              << "EVALUATING RES: PERTURBED" << std::endl
              << "--------------------------------------------------" << std::endl;
    pk->FunctionalResidual(0., 8640, soln_old, soln, res2);
    db->WriteVector("mat_diff rhs", pk_richards->get_operator()->rhs().ptr(), true);
    pk->ErrorNorm(soln, res2);

    // Compute r(p + dp) - r(p) ~= J * dp
    // -- diff residuals
    res2->Update(-1, *res, 1.);

    // -- apply the PC
    std::cout << std::endl
              << "APPLY PC" << std::endl
              << "--------------------------------------------------" << std::endl;
    res->PutScalar(0.);
    pk->preconditioner()->Apply(*dsoln->SubVector(0)->Data(), *res->SubVector(0)->Data());

    db->WriteVector("dres", res2->SubVector(0)->Data().ptr(), true);
    db->WriteVector("Jdp", res->SubVector(0)->Data().ptr(), true);

    // -- subtract J dp
    double norm;
    res->Update(-1 / eps, *res2, 1. / eps);
    res->NormInf(&norm);
    CHECK_CLOSE(0., norm, 1.e-3);
    if (norm > 1.e-3) {
      std::cout << "res(p + dp) - res(p) - J(dp):" << std::endl;
      res->Print(std::cout);
    }
  };
}
