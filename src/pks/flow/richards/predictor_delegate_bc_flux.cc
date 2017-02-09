/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Delegate for modifying the predictor in the case of infiltration into dry soil.

  License: BSD
  Authors: Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include <boost/math/tools/roots.hpp>
#include "predictor_delegate_bc_flux.hh"

#include "Op.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1

bool PredictorDelegateBCFlux::ModifyPredictor(const Teuchos::Ptr<CompositeVector>& u) {
  Epetra_MultiVector& u_f = *u->ViewComponent("face",false);

  int nfaces = bc_values_->size();
  for (int f=0; f!=nfaces; ++f) {
    if ((*bc_markers_)[f] == Operators::OPERATOR_BC_NEUMANN) {
      double lambda = u_f[0][f];
      // only do if below saturated
      if (lambda < 101325.) {
        int ierr = CalculateLambdaToms_(f, u, lambda);
        ASSERT(!ierr);
        if (!ierr) u_f[0][f] = lambda;
      }
    }
  }
  return true;
}


Teuchos::RCP<PredictorDelegateBCFlux::FluxBCFunctor>
PredictorDelegateBCFlux::CreateFunctor_(int f,
        const Teuchos::Ptr<const CompositeVector>& pres) {
  // inner cell
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  ASSERT(cells.size() == 1);
  int c = cells[0];

  // that cell's faces
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;
  mesh_->cell_get_faces_and_dirs(c, &faces, &dirs);

  // index within that cell's faces
  unsigned int n = std::find(faces.begin(), faces.end(), f) - faces.begin();
  ASSERT(n != faces.size());

  Teuchos::RCP< std::vector<double> > Aff = Teuchos::rcp(new std::vector<double>());
  Aff->resize(faces.size());
  Teuchos::RCP< std::vector<double> > lambda = Teuchos::rcp(new std::vector<double>());
  lambda->resize(faces.size());

  // unscale the Aff for my cell with the rel perm that was used to calculate it
  double Krel = (*S_next_->GetFieldData("upwind_relative_permeability")
                 ->ViewComponent("face",false))[0][f];

  // fill the arrays
  for (unsigned int i=0; i!=faces.size(); ++i) {
    (*Aff)[i] = matrix_->local_matrices()->matrices[c](n,i) / Krel;
    (*lambda)[i] = (*pres)("face",faces[i]);
  }

  // gravity flux
  double bc_flux = mesh_->face_area(f) * (*bc_values_)[f];
  double gflux = ((*matrix_->global_operator()->rhs()->ViewComponent("face",false))[0][faces[n]] + bc_flux) / Krel;

#if DEBUG_FLAG
  std::cout << "   Aff = ";
  for (unsigned int i=0; i!=faces.size(); ++i) std::cout << (*Aff)[i] << ", ";
  std::cout << std::endl << "   lambda = ";
  for (unsigned int i=0; i!=faces.size(); ++i) std::cout << (*lambda)[i] << ", ";
  std::cout << std::endl << "   p_cell = " << (*pres)("cell",c) << std::endl;
  std::cout << "    and init K_rel = " << wrms_->second[(*wrms_->first)[c]]->k_relative(wrms_->second[(*wrms_->first)[c]]->saturation(101325. - (*lambda)[n])) << std::endl;
  std::cout << "    to match fluxes: bc = " << bc_flux << " and grav = " << gflux << std::endl;
#endif

  // create and return
  return Teuchos::rcp(new FluxBCFunctor(Aff, lambda, n, (*pres)("cell",c),
          bc_flux, gflux, dirs[n], 101325.0, wrms_->second[(*wrms_->first)[c]]));
}

int PredictorDelegateBCFlux::CalculateLambdaToms_(int f,
        const Teuchos::Ptr<const CompositeVector>& pres, double& lambda) {

#if DEBUG_FLAG
  std::cout << " Flux correcting face " << f << ": q = " << (*bc_values_)[f] << std::endl;
#endif

  // start by making sure lambda is a reasonable guess, which may not be the case
  if (std::abs(lambda) > 1.e7) lambda = 101325.;

  Teuchos::RCP<FluxBCFunctor> func = CreateFunctor_(f,pres);

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

  if (res == 0) return 0;
    

  if (res > 0.) {
    left = lambda;
    lres = res;
    right = std::max(lambda,101325.);
    rres = (*func)(right);
    while (rres > 0.) {
      right += 101325.;
      rres = (*func)(right);
    }

  } else {
    right = lambda;
    rres = res;

    std::cout << "RIGHT = " << right << ", " << rres << std::endl;
    left = std::min(101325., lambda);
    lres = (*func)(left);
    while (lres < 0.) {
      std::cout << "LEFT = " << left << ", " << lres << std::endl;
      left -= 101325.;
      lres = (*func)(left);
    }
  }

#if DEBUG_FLAG
  std::cout << "   bracket (res): " << left << " (" << lres << "), "
            << right << " (" << rres << ")" << std::endl;
#endif

  std::pair<double,double> result =
      boost::math::tools::toms748_solve(*func, left, right, lres, rres, tol, actual_it);
  if (actual_it >= max_it) {
    std::cout << " Failed to converged in " << actual_it << " steps." << std::endl;
    return 3;
  }

  lambda = (result.first + result.second) / 2.;
#if DEBUG_FLAG
  std::cout << "  Converged to " << lambda << " in " << actual_it << " steps." << std::endl;

  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, AmanziMesh::USED, &cells);
  ASSERT(cells.size() == 1);
  int c = cells[0];

  std::cout << "      with k_rel = " << wrms_->second[(*wrms_->first)[c]]->k_relative(wrms_->second[(*wrms_->first)[c]]->saturation(101325. - lambda)) << std::endl;


#endif
  return 0;
}



} // namespace
} // namespace
