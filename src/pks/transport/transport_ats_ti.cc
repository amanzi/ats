/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Transport PK

*/

#include <algorithm>

#include "OperatorDefs.hh"
#include "transport_ats.hh"

namespace Amanzi {
namespace Transport {

/* *******************************************************************
 * Routine takes a parallel overlapping vector C and returns a parallel
 * overlapping vector F(C).
 ****************************************************************** */
void
Transport_ATS::FunctionalTimeDerivative(double t,
                                        const Epetra_Vector& component,
                                        Epetra_Vector& f_component)
{
  int nfaces_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  // distribute vector
  auto component_tmp = Teuchos::rcp(new Epetra_Vector(component));
  component_tmp->Import(component, tcc->importer("cell"), Insert);

  Teuchos::ParameterList recon_list = plist_->sublist("reconstruction");
  lifting_->Init(recon_list);
  lifting_->Compute(component_tmp);

  // extract boundary conditions for the current component
  std::vector<int> bc_model(nfaces_all, Operators::OPERATOR_BC_NONE);
  std::vector<double> bc_value(nfaces_all);

  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;
          std::vector<double>& values = it->second;

          bc_model[f] = Operators::OPERATOR_BC_DIRICHLET;
          bc_value[f] = values[i];
        }
      }
    }
  }

  Teuchos::RCP<const Epetra_MultiVector> flux =
    S_->Get<CompositeVector>(flux_key_, Tags::NEXT).ViewComponent("face", true);
  limiter_->Init(recon_list, flux);
  limiter_->ApplyLimiter(component_tmp, 0, lifting_, bc_model, bc_value);
  lifting_->data()->ScatterMasterToGhosted("cell");

  // ADVECTIVE FLUXES
  // We assume that limiters made their job up to round-off errors.
  // Min-max condition will enforce robustness w.r.t. these errors.

  f_component.PutScalar(0.0);
  for (int f = 0; f < nfaces_all; f++) { // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    double u1, u2, umin, umax;
    if (c1 >= 0 && c2 >= 0) {
      u1 = (*component_tmp)[c1];
      u2 = (*component_tmp)[c2];
      umin = std::min(u1, u2);
      umax = std::max(u1, u2);
    } else if (c1 >= 0) {
      u1 = u2 = umin = umax = (*component_tmp)[c1];
    } else if (c2 >= 0) {
      u1 = u2 = umin = umax = (*component_tmp)[c2];
    }

    double u = std::abs((*flux)[0][f]);
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);

    double upwind_tcc, tcc_flux;
    if (c1 >= 0 && c1 < f_component.MyLength() && c2 >= 0 && c2 < f_component.MyLength()) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;
      f_component[c2] += tcc_flux;

    } else if (c1 >= 0 && c1 < f_component.MyLength() && (c2 >= f_component.MyLength() || c2 < 0)) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;

    } else if (c1 >= 0 && c1 < f_component.MyLength() && (c2 < 0)) {
      upwind_tcc = component[c1];
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c1] -= tcc_flux;

    } else if (c1 >= f_component.MyLength() && c2 >= 0 && c2 < f_component.MyLength()) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = u * upwind_tcc;
      f_component[c2] += tcc_flux;
    }
  }

  // process external sources
  if (srcs_.size() != 0) {
    ComputeAddSourceTerms_(t, 1., f_component, current_component_, current_component_);
  }

  S_->GetEvaluator(porosity_key_, Tags::NEXT).Update(*S_, name_);
  const Epetra_MultiVector& phi =
    *S_->Get<CompositeVector>(porosity_key_, Tags::NEXT).ViewComponent("cell", false);
  for (int c = 0; c < f_component.MyLength(); c++) { // calculate conservative quantatity
    double vol_phi_ws_den =
      mesh_->getCellVolume(c) * phi[0][c] * (*ws_prev_)[0][c] * (*mol_dens_prev_)[0][c];
    if ((*ws_prev_)[0][c] < 1e-12)
      vol_phi_ws_den = mesh_->getCellVolume(c) * phi[0][c] * (*ws_)[0][c] * (*mol_dens_)[0][c];

    if (vol_phi_ws_den > water_tolerance_) { f_component[c] /= vol_phi_ws_den; }
  }

  // boundary conditions for advection
  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (current_component_ == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;
          std::vector<double>& values = it->second;
          int c2 = (*downwind_cell_)[f];

          if (c2 >= 0 && f < nfaces_owned) {
            double u = std::abs((*flux)[0][f]);
            double vol_phi_ws_den =
              mesh_->getCellVolume(c2) * phi[0][c2] * (*ws_prev_)[0][c2] * (*mol_dens_prev_)[0][c2];
            if ((*ws_prev_)[0][c2] < 1e-12)
              vol_phi_ws_den =
                mesh_->getCellVolume(c2) * phi[0][c2] * (*ws_)[0][c2] * (*mol_dens_)[0][c2];

            double tcc_flux = u * values[i];
            if (vol_phi_ws_den > water_tolerance_) { f_component[c2] += tcc_flux / vol_phi_ws_den; }
          }
        }
      }
    }
  }
}


} // namespace Transport
} // namespace Amanzi
