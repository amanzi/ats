/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

/*
  Delegate for modifying the predictor in the case of infiltration into dry soil.

*/

#include "predictor_delegate_bc_flux.hh"

#include "Op.hh"
#include "Brent.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 0

bool
PredictorDelegateBCFlux::ModifyPredictor(const Teuchos::Ptr<CompositeVector>& u)
{
  auto u_f = u->viewComponent<MemSpace_kind::HOST>("face", false);
  auto markers_f = bc_markers_->viewComponent<MemSpace_kind::HOST>("face", false);

  int nfaces = u_f.extent(0);
  for (auto& region_wrm : wrms_) {
    auto rfaces = mesh_->getSetEntities<MemSpace_kind::HOST>(
      region_wrm.first, AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);

    for (AmanziMesh::Entity_ID f : rfaces) {
      if (markers_f(f, 0) == Operators::OPERATOR_BC_NEUMANN) {
        double lambda = u_f(f, 0);
        // only do if below saturated
        if (lambda < 101325.) {
          auto fcells = mesh_->getFaceCells(f);
          AMANZI_ASSERT(fcells.size() == 1);
          int ierr = CalculateLambda_(f, fcells(0), region_wrm.second, u, lambda);
          AMANZI_ASSERT(!ierr);
          if (!ierr) u_f(f, 0) = lambda;
        }
      }
    }
  }
  return true;
}


Teuchos::RCP<PredictorDelegateBCFlux::FluxBCFunctor>
PredictorDelegateBCFlux::CreateFunctor_(AmanziMesh::Entity_ID f,
                                        AmanziMesh::Entity_ID c,
                                        const Teuchos::RCP<WRMModel_type>& wrm,
                                        const Teuchos::Ptr<const CompositeVector>& pres)
{
  // get cell's faces
  auto [faces, dirs] = mesh_->getCellFacesAndDirections(c);

  // index within that cell's faces
  unsigned int n = 0;
  for (; n != faces.size(); ++n)
    if (faces(n) == f) break;
  AMANZI_ASSERT(n != faces.size());

  // local matrix row, vector
  Kokkos::View<double*, Amanzi::DefaultHost> Aff("Aff", faces.size());
  Kokkos::View<double*, Amanzi::DefaultHost> lambda("lambda", faces.size());

  // collect physics
  const auto& pres_f = pres->viewComponent<MemSpace_kind::HOST>("face", false);
  const auto& pres_c = pres->viewComponent<MemSpace_kind::HOST>("cell", false);
  const auto& rhs_f = matrix_->global_operator()->rhs()->viewComponent<MemSpace_kind::HOST>("face", false);
  const auto& bc_values = bc_values_->viewComponent<MemSpace_kind::HOST>("face", false);

  // unscale the Aff for my cell with rel perm
  double Krel = wrm->getModel().k_relative(wrm->getModel().saturation(101325. - pres_f(f, 0)));

  // fill the arrays
  const auto Aff_g = matrix_->local_op()->A.at_host(c);
  for (unsigned int i = 0; i != faces.size(); ++i) {
    Aff(i) = Aff_g(n, i) / Krel;
    lambda(i) = pres_f(faces(i), 0);
  }

  // gravity flux
  double bc_flux = mesh_->getFaceArea(f) * bc_values(f, 0);
  double gflux = rhs_f(faces(n), 0) / Krel;

  // #if DEBUG_FLAG
  //   std::cout << "   Aff = ";
  //   for (unsigned int i = 0; i != faces.size(); ++i) std::cout << (*Aff)[i] << ", ";
  //   std::cout << std::endl << "   lambda = ";
  //   for (unsigned int i = 0; i != faces.size(); ++i) std::cout << (*lambda)[i] << ", ";
  //   std::cout << std::endl << "   p_cell = " << (*pres)("cell", c) << std::endl;
  //   std::cout << "    and init K_rel = "
  //             << wrms_->second[(*wrms_->first)[c]]->k_relative(
  //                  wrms_->second[(*wrms_->first)[c]]->saturation(101325. - (*lambda)[n]))
  //             << std::endl;
  //   std::cout << "    to match fluxes: bc = " << bc_flux << " and grav = " << gflux << std::endl;
  // #endif

  // create and return
  return Teuchos::rcp(new FluxBCFunctor(
    Aff, lambda, n, pres_c(c, 0), bc_flux, gflux, dirs(n), 101325.0, wrm->getModel()));
}

int
PredictorDelegateBCFlux::CalculateLambda_(AmanziMesh::Entity_ID f,
                                          AmanziMesh::Entity_ID c,
                                          const Teuchos::RCP<WRMModel_type>& wrm,
                                          const Teuchos::Ptr<const CompositeVector>& pres,
                                          double& lambda)
{
  // #if DEBUG_FLAG
  //   std::cout << " Flux correcting face " << f << ": q = " << (*bc_values_)[f] << std::endl;
  // #endif

  // start by making sure lambda is a reasonable guess, which may not be the case
  if (std::abs(lambda) > 1.e7) lambda = 101325.;

  Teuchos::RCP<FluxBCFunctor> func = CreateFunctor_(f, c, wrm, pres);

  // -- convergence criteria
  const auto& bc_values = bc_values_->viewComponent<MemSpace_kind::HOST>("face", false);
  double eps = std::max(1.e-4 * std::abs(bc_values(f, 0)), 1.e-8);
  int max_it = 100;
  int actual_it(max_it);

  double res = (*func)(lambda);
  double left = 0.;
  double right = 0.;
  double lres = 0.;
  double rres = 0.;

  if (std::abs(res) < eps) {
#if DEBUG_FLAG
    std::cout << "  Converged to " << lambda << " in " << 0 << " steps." << std::endl;
#endif
    return 0;
  }


  if (res > 0.) {
    left = lambda;
    lres = res;
    right = std::max(lambda, 101325.);
    rres = (*func)(right);
    while (rres > 0.) {
      right += 101325.;
      rres = (*func)(right);
    }

  } else {
    right = lambda;
    rres = res;

    left = std::min(101325., lambda);
    lres = (*func)(left);
    while (lres < 0.) {
      left -= 101325.;
      lres = (*func)(left);
    }
  }

#if DEBUG_FLAG
  std::cout << "   bracket (res): " << left << " (" << lres << "), " << right << " (" << rres << ")"
            << std::endl;
#endif

  lambda = Utils::findRootBrent(*func, left, right, eps, &actual_it);
  if (actual_it >= max_it) {
    std::cout << " Failed to converged in " << actual_it << " steps." << std::endl;
    return 3;
  }

#if DEBUG_FLAG
  std::cout << "  Converged to " << lambda << " in " << actual_it << " steps." << std::endl;

  auto cells = mesh_->getFaceCells(f);
  AMANZI_ASSERT(cells.size() == 1);
  int c = cells[0];

  std::cout << "      with k_rel = "
            << wrms_->second[(*wrms_->first)[c]]->k_relative(
                 wrms_->second[(*wrms_->first)[c]]->saturation(101325. - lambda))
            << std::endl;


#endif
  return 0;
}


} // namespace Flow
} // namespace Amanzi
