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

#include "boost/math/tools/roots.hpp"
#include "predictor_delegate_bc_flux.hh"

#include "Op.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 0

bool
PredictorDelegateBCFlux::ModifyPredictor(const Teuchos::Ptr<CompositeVector>& u)
{
  Epetra_MultiVector& u_f = *u->ViewComponent("face", false);

  int nfaces = bc_values_->size();
  for (int f = 0; f != nfaces; ++f) {
    if ((*bc_markers_)[f] == Operators::OPERATOR_BC_NEUMANN) {
      double lambda = u_f[0][f];
      // only do if below saturated
      if (lambda < 101325.) {
        int ierr = CalculateLambdaToms_(f, u, lambda);
        AMANZI_ASSERT(!ierr);
        if (!ierr) u_f[0][f] = lambda;
      }
    }
  }
  return true;
}


Teuchos::RCP<PredictorDelegateBCFlux::FluxBCFunctor>
PredictorDelegateBCFlux::CreateFunctor_(int f, const Teuchos::Ptr<const CompositeVector>& pres)
{
  // inner cell and its water retention model
  auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
  AMANZI_ASSERT(cells.size() == 1);
  int c = cells[0];

  // that cell's faces
  const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);

  // index within that cell's faces
  unsigned int n = std::find(faces.begin(), faces.end(), f) - faces.begin();
  AMANZI_ASSERT(n != faces.size());

  // local matrix row, vector
  Teuchos::RCP<std::vector<double>> Aff = Teuchos::rcp(new std::vector<double>());
  Aff->resize(faces.size());
  Teuchos::RCP<std::vector<double>> lambda = Teuchos::rcp(new std::vector<double>());
  lambda->resize(faces.size());

  // collect physics
  const auto& wrm = wrms_->second[(*wrms_->first)[c]];
  const auto& pres_f = *pres->ViewComponent("face", false);
  const auto& pres_c = *pres->ViewComponent("cell", false);
  const auto& rhs_f = *matrix_->global_operator()->rhs()->ViewComponent("face", false);

  // unscale the Aff for my cell with rel perm
  double Krel = wrm->k_relative(wrm->saturation(101325. - pres_f[0][f]));

  // fill the arrays
  const auto& Aff_g = matrix_->local_op()->matrices[c];
  for (unsigned int i = 0; i != faces.size(); ++i) {
    (*Aff)[i] = Aff_g(n, i) / Krel;
    (*lambda)[i] = pres_f[0][faces[i]];
  }

  // gravity flux
  double bc_flux = mesh_->getFaceArea(f) * (*bc_values_)[f];
  double gflux = rhs_f[0][faces[n]] / Krel;

#if DEBUG_FLAG
  std::cout << "   Aff = ";
  for (unsigned int i = 0; i != faces.size(); ++i) std::cout << (*Aff)[i] << ", ";
  std::cout << std::endl << "   lambda = ";
  for (unsigned int i = 0; i != faces.size(); ++i) std::cout << (*lambda)[i] << ", ";
  std::cout << std::endl << "   p_cell = " << (*pres)("cell", c) << std::endl;
  std::cout << "    and init K_rel = "
            << wrms_->second[(*wrms_->first)[c]]->k_relative(
                 wrms_->second[(*wrms_->first)[c]]->saturation(101325. - (*lambda)[n]))
            << std::endl;
  std::cout << "    to match fluxes: bc = " << bc_flux << " and grav = " << gflux << std::endl;
#endif

  // create and return
  return Teuchos::rcp(
    new FluxBCFunctor(Aff, lambda, n, pres_c[0][c], bc_flux, gflux, dirs[n], 101325.0, wrm));
}

int
PredictorDelegateBCFlux::CalculateLambdaToms_(int f,
                                              const Teuchos::Ptr<const CompositeVector>& pres,
                                              double& lambda)
{
#if DEBUG_FLAG
  std::cout << " Flux correcting face " << f << ": q = " << (*bc_values_)[f] << std::endl;
#endif

  // start by making sure lambda is a reasonable guess, which may not be the case
  if (std::abs(lambda) > 1.e7) lambda = 101325.;

  Teuchos::RCP<FluxBCFunctor> func = CreateFunctor_(f, pres);

  // -- convergence criteria
  double eps = std::max(1.e-4 * std::abs((*bc_values_)[f]), 1.e-8);
  Tol_ tol(eps);
  boost::uintmax_t max_it = 100;
  boost::uintmax_t actual_it(max_it);

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

  std::pair<double, double> result =
    boost::math::tools::toms748_solve(*func, left, right, lres, rres, tol, actual_it);
  if (actual_it >= max_it) {
    std::cout << " Failed to converged in " << actual_it << " steps." << std::endl;
    return 3;
  }

  lambda = (result.first + result.second) / 2.;
#if DEBUG_FLAG
  std::cout << "  Converged to " << lambda << " in " << actual_it << " steps." << std::endl;

  AmanziMesh::Entity_ID_List cells;
  cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
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
