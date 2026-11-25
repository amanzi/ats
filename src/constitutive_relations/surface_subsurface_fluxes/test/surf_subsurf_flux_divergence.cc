/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

/*
  FluxDivergenceEvaluator must give the same cell divergence regardless of how
  the mesh is partitioned.
*/

#include <cmath>
#include <iostream>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "EvaluatorPrimary.hh"
#include "State.hh"

#include "flux_divergence_evaluator.hh"

using namespace Amanzi;

namespace {

// specify nonzero and spatially varying divergence
AmanziGeometry::Point
velocity(const AmanziGeometry::Point& x)
{
  return AmanziGeometry::Point(1.0 + x[0] * x[0], 0.5 + 0.25 * x[1]);
}

void
evaluateDivergence(const Teuchos::RCP<AmanziMesh::Mesh>& mesh, State& S)
{
  S.RegisterDomainMesh(mesh);

  // face fluxes
  Teuchos::ParameterList flux_list("water_flux");
  flux_list.set<std::string>("tag", "");
  S.Require<CompositeVector, CompositeVectorSpace>("water_flux", Tags::DEFAULT, "water_flux")
    .SetMesh(mesh)
    ->SetComponent("face", AmanziMesh::Entity_kind::FACE, 1);
  auto flux_eval =
    Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(flux_list));
  S.SetEvaluator("water_flux", Tags::DEFAULT, flux_eval);

  // flux divergence
  Teuchos::ParameterList div_list("water_flux_divergence");
  div_list.set<std::string>("tag", "");
  S.Require<CompositeVector, CompositeVectorSpace>("water_flux_divergence", Tags::DEFAULT)
    .SetMesh(mesh)
    ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  auto div_eval = Teuchos::rcp(new Relations::FluxDivergenceEvaluator(div_list));
  S.SetEvaluator("water_flux_divergence", Tags::DEFAULT, div_eval);

  S.Setup();
  S.GetW<CompositeVector>("water_flux", Tags::DEFAULT, "water_flux").PutScalarMasterAndGhosted(0.);
  S.GetRecordW("water_flux", Tags::DEFAULT, "water_flux").set_initialized();
  S.Initialize();

  auto& flux = *S.GetW<CompositeVector>("water_flux", Tags::DEFAULT, "water_flux")
                  .ViewComponent("face", false);
  for (int f = 0; f != flux.MyLength(); ++f) {
    flux[0][f] = velocity(mesh->getFaceCentroid(f)) * mesh->getFaceNormal(f);
  }
  flux_eval->SetChanged();

  S.GetEvaluator("water_flux_divergence", Tags::DEFAULT).Update(S, "test");
}

} // namespace


SUITE(SURF_SUBSURF_FLUX_DIVERGENCE)
{
  TEST(PARTITION_INDEPENDENCE)
  {
    auto comm = getDefaultComm();
    auto comm_self = getCommSelf();

    // 8 x 2 surface-like 2D mesh: on 2 ranks the partition boundary cuts
    // through faces that carry nonzero flux
    const double x1 = 8.0, y1 = 2.0;
    const int nx = 8, ny = 2;

    State S_par, S_ser;
    auto mesh_par = AmanziMesh::MeshFactory(comm).create(0.0, 0.0, x1, y1, nx, ny);
    auto mesh_ser = AmanziMesh::MeshFactory(comm_self).create(0.0, 0.0, x1, y1, nx, ny);
    evaluateDivergence(mesh_par, S_par);
    evaluateDivergence(mesh_ser, S_ser);

    const auto& div_par = *S_par.Get<CompositeVector>("water_flux_divergence", Tags::DEFAULT)
                             .ViewComponent("cell", false);
    const auto& div_ser = *S_ser.Get<CompositeVector>("water_flux_divergence", Tags::DEFAULT)
                             .ViewComponent("cell", false);
    CHECK_EQUAL(nx * ny, div_ser.MyLength());

    // every owned cell of the distributed mesh must match the equivalent serial cell
    int n_mismatch = 0;
    for (int c = 0; c != div_par.MyLength(); ++c) {
      const auto& xc = mesh_par->getCellCentroid(c);
      int cs = -1;
      for (int s = 0; s != div_ser.MyLength(); ++s) {
        if (AmanziGeometry::norm(mesh_ser->getCellCentroid(s) - xc) < 1.e-10) {
          cs = s;
          break;
        }
      }
      CHECK(cs >= 0);
      if (cs < 0) continue;

      double scale = std::max(1.0, std::abs(div_ser[0][cs]));
      if (std::abs(div_par[0][c] - div_ser[0][cs]) > 1.e-12 * scale) {
        ++n_mismatch;
        std::cout << "rank " << comm->MyPID() << ": cell at (" << xc[0] << ", " << xc[1]
                  << ") divergence " << div_par[0][c] << " (distributed) != " << div_ser[0][cs]
                  << " (serial)" << std::endl;
      }
    }

    int n_mismatch_global = 0;
    comm->SumAll(&n_mismatch, &n_mismatch_global, 1);
    CHECK_EQUAL(0, n_mismatch_global);

    // sanity check - ensure div is not zero
    double max_div = 0.;
    for (int s = 0; s != div_ser.MyLength(); ++s) max_div = std::max(max_div, std::abs(div_ser[0][s]));
    CHECK(max_div > 1.e-3);
  }
}
