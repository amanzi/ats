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

//
// Given C at t_old, computes cq_flux += div q * C * dt
//
// NOTE, we can assume safely that most normal things (water contents, flux,
// etc) are already updated, and upwind/downwind cells are already identified,
// because we have called ComputeStableTimestep already.
//
// This uses a first-order donor upwind scheme.
//
// DEV NOTE: this requires that tcc is ghosted and scattered, and that cq_flux is
// owned.
void
Transport_ATS::AddAdvection_FirstOrderUpwind_(double t_old,
        double t_new,
        const Epetra_MultiVector& tcc,
        Epetra_MultiVector& cq_flux)
{
  double dt = t_new - t_old;

  int nfaces_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);
  const Epetra_MultiVector& flux = *S_->Get<CompositeVector>(flux_key_, tag_next_)
    .ViewComponent("face", true);

  // advance all components at once
  for (int f = 0; f < nfaces_all; f++) {
    // flow moves from upwind cell (c1) to downwind cell (c2).
    // If ci < 0 || ci > ncells_owned -> indicates boundary or halo cells (i=1,2)
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];
    double u = std::abs(flux[0][f]);

    if (c1 >= 0 && c1 < cq_flux.MyLength() && c2 >= 0 && c2 < cq_flux.MyLength()) {
      // Here c1 & c2 are inside local domain. Update solute fluxes for both cells
      for (int i = 0; i < num_aqueous_; i++) {
        double tcc_flux = dt * u * tcc[i][c1];
        cq_flux[i][c1] -= tcc_flux;
        cq_flux[i][c2] += tcc_flux;
      }
      // Update (tag next) water fluxes for both cells
      cq_flux[num_aqueous_ + 1][c1] -= dt * u;
      cq_flux[num_aqueous_ + 1][c2] += dt * u;

    } else if (c1 >= 0 && c1 < cq_flux.MyLength() && (c2 >= cq_flux.MyLength() || c2 < 0)) {
      // downind cell c2 is boundary or belong to another domain owned by other processors
      // Update solute flux for c1
      for (int i = 0; i < num_aqueous_; i++) {
        double tcc_flux = dt * u * tcc[i][c1];
        cq_flux[i][c1] -= tcc_flux;
      }
      // Update or subtract (tag next) water fluxes for c1
      cq_flux[num_aqueous_ + 1][c1] -= dt * u;

    } else if (c1 >= cq_flux.MyLength() && c2 >= 0 && c2 < cq_flux.MyLength()) {
      // upwind cell c1 is boundary or belong to another domain owned by other processors
      // Update solute flux for c2
      for (int i = 0; i < num_aqueous_; i++) {
        double tcc_flux = dt * u * tcc[i][c1];
        cq_flux[i][c2] += tcc_flux;
      }
      // Update or add (tag next) water fluxes for c2
      cq_flux[num_aqueous_ + 1][c2] += dt * u;

    } else if (c2 < 0 && c1 >= 0 && c1 < cq_flux.MyLength()) {
      // this case is very similar to line 1165, except c2<0 only. Why we don't update solute flux??? --PL
      //
      // Negative cell value implies the face is a domain boundary (not process
      // boundary).  Solute fluxes over domain boundaries (to c1) are taken
      // care of in the below boundary loop, but water fluxes aren't because
      // water fluxes may not bring solute, so deal with water fluxes
      // here. --ETC
      cq_flux[num_aqueous_ + 1][c1] -= dt * u;

    } else if (c1 < 0 && c2 >= 0 && c2 < cq_flux.MyLength()) {
      // why no solute update??? --PL
      cq_flux[num_aqueous_ + 1][c2] += dt * u;
    }
  }

  // loop over exterior boundary sets
  //
  // Why no check on the boundary type?  This can be Dirichlet or Neumann?  Or
  // is this hard-coded as just Dirichlet data? --ETC
  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
      int f = it->first;

      std::vector<double>& values = it->second;
      int c2 = (*downwind_cell_)[f];
      int c1 = (*upwind_cell_)[f];

      double u = std::abs(flux[0][f]);
      if (c2 >= 0 && c2 < cq_flux.MyLength()) {
        for (int i = 0; i < ncomp; i++) {
          int k = tcc_index[i];
          if (k < num_aqueous_) {
            double tcc_flux = dt * u * values[i];
            cq_flux[k][c2] += tcc_flux;
          }
        }
      }
    }
  }
}


void
Transport_ATS::AddAdvection_SecondOrderUpwind_(double t_old,
        double t_new,
        const Epetra_MultiVector& tcc,
        Epetra_MultiVector& cq_flux)
{
  for (int i = 0; i != num_aqueous_; ++i) {
    AddAdvection_SecondOrderUpwind_(t_old, t_new, *tcc(i), *cq_flux(i), i);
  }
}


//
// Given C_old, this computes div Cq using a second order limiter
//
// DEV NOTE: requires that tcc is ghosted and scattered, while cq_flux is
// owned.
void
Transport_ATS::AddAdvection_SecondOrderUpwind_(double t_old,
        double t_new,
        const Epetra_Vector& tcc,
        Epetra_Vector& cq_flux,
        int component)
{
  double dt = t_new - t_old;

  int nfaces_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::OWNED);
  int nfaces_all =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::FACE, AmanziMesh::Parallel_kind::ALL);

  // API funkiness
  Teuchos::RCP<const Epetra_MultiVector> tcc_rcp = Teuchos::rcpFromRef<const Epetra_MultiVector>(tcc);
  lifting_->Compute(tcc_rcp, 0);

  // populate boundary conditions for the current component
  PopulateBoundaryData_(component, *adv_bcs_);
  auto& bc_model = adv_bcs_->bc_model();
  auto& bc_value = adv_bcs_->bc_value();

  auto flux = S_->Get<CompositeVector>(flux_key_, tag_next_).ViewComponent("face", true);
  limiter_->SetFlux(flux);
  limiter_->ApplyLimiter(tcc_rcp, 0, lifting_, bc_model, bc_value);

  // ADVECTIVE FLUXES
  // We assume that limiters made their job up to round-off errors.
  // Min-max condition will enforce robustness w.r.t. these errors.
  for (int f = 0; f < nfaces_all; f++) { // loop over master and slave faces
    int c1 = (*upwind_cell_)[f];
    int c2 = (*downwind_cell_)[f];

    double u1, u2, umin, umax;
    if (c1 >= 0 && c2 >= 0) {
      u1 = tcc[c1];
      u2 = tcc[c2];
      umin = std::min(u1, u2);
      umax = std::max(u1, u2);
    } else if (c1 >= 0) {
      u1 = u2 = umin = umax = tcc[c1];
    } else if (c2 >= 0) {
      u1 = u2 = umin = umax = tcc[c2];
    }

    double u = std::abs((*flux)[0][f]);
    const AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);

    double upwind_tcc, tcc_flux;
    if (c1 >= 0 && c1 < cq_flux.MyLength() && c2 >= 0 && c2 < cq_flux.MyLength() ) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = dt * u * upwind_tcc;
      cq_flux[c1] -= tcc_flux;
      cq_flux[c2] += tcc_flux;

    } else if (c1 >= 0 && c1 < cq_flux.MyLength() && (c2 >= cq_flux.MyLength() || c2 < 0)) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = dt * u * upwind_tcc;
      cq_flux[c1] -= tcc_flux;

    } else if (c1 >= 0 && c1 < cq_flux.MyLength() && (c2 < 0)) {
      upwind_tcc = tcc[c1];
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = dt * u * upwind_tcc;
      cq_flux[c1] -= tcc_flux;

    } else if (c1 >= cq_flux.MyLength() && c2 >= 0 && c2 < cq_flux.MyLength() ) {
      upwind_tcc = limiter_->getValue(c1, xf);
      upwind_tcc = std::max(upwind_tcc, umin);
      upwind_tcc = std::min(upwind_tcc, umax);

      tcc_flux = dt * u * upwind_tcc;
      cq_flux[c2] += tcc_flux;
    }
  }

  // boundary conditions for advection
  for (int m = 0; m < bcs_.size(); m++) {
    std::vector<int>& tcc_index = bcs_[m]->tcc_index();
    int ncomp = tcc_index.size();

    for (int i = 0; i < ncomp; i++) {
      if (component == tcc_index[i]) {
        for (auto it = bcs_[m]->begin(); it != bcs_[m]->end(); ++it) {
          int f = it->first;
          std::vector<double>& values = it->second;
          int c2 = (*downwind_cell_)[f];

          if (c2 >= 0 && f < nfaces_owned) {
            double u = std::abs((*flux)[0][f]);
            double tcc_flux = dt * u * values[i];
            cq_flux[c2] += tcc_flux;
          }
        }
      }
    }
  }
}



// Given the conserved quantity, computes new tcc
//
// DEV NOTE: all vectors are owned only
//
void
Transport_ATS::InvertTccNew_(const Epetra_MultiVector& conserve_qty,
                             Epetra_MultiVector& tcc,
                             Epetra_MultiVector* solid_qty,
                             bool include_current_water_mass)
{
  // note this was updated in StableStep
  const Epetra_MultiVector& lwc_new = *S_->Get<CompositeVector>(lwc_key_, tag_next_)
    .ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S_->Get<CompositeVector>(cv_key_, tag_next_)
    .ViewComponent("cell", false);

  // calculate the new conc based on advected term
  for (int c = 0; c < conserve_qty.MyLength(); c++) {
    double water_new = lwc_new[0][c];
    double water_sink = conserve_qty[num_aqueous_][c];
    double water_total = water_sink + water_new;

    for (int i = 0; i != num_components_; ++i) {
      double tcc_max = tcc_max_[component_names_[i]];

      if (water_new > water_tolerance_ / cv[0][c] && conserve_qty[i][c] > 0) {
        // there is both water and stuff present at the new time this is stuff
        // at the new time + stuff leaving through the domain coupling, divided
        // by water of both
        tcc[i][c] = conserve_qty[i][c] / water_total;
        if (include_current_water_mass) conserve_qty[num_aqueous_][c] = water_total;

        if (solid_qty && tcc_max > 0 && tcc[i][c] > tcc_max) {
          (*solid_qty)[i][c] += (tcc[i][c] - tcc_max) * water_total;
          tcc[i][c] = tcc_max;
        }

      } else if (water_sink > water_tolerance_ / cv[0][c] && conserve_qty[i][c] > 0) {
        // there is water and stuff leaving through the domain coupling, but it
        // all leaves (none at the new time)
        tcc[i][c] = 0.;
      } else {
        // there is no water leaving, and no water at the new time.  Change any
        // stuff into solid
        if (solid_qty) (*solid_qty)[i][c] += std::max(conserve_qty[i][c], 0.);
        conserve_qty[i][c] = 0.;
        tcc[i][c] = 0.;
      }

      // dissolve solid
      if (solid_qty &&
          water_total > water_tolerance_ / cv[0][c] && // water available for dissolution
          (*solid_qty)[i][c] > 0. && // solid available to dissolve
          tcc[i][c] < tcc_max) {
        double dissolved_qty = std::min( (tcc_max - tcc[i][c]) * water_total, (*solid_qty)[i][c]);
        if (water_new > water_tolerance_ / cv[0][c]) {
          // new water -- dissolve in place
          tcc[i][c] += dissolved_qty / water_new;
        } // otherwise dissolve into water_sink, advect out, tcc stays 0

        (*solid_qty)[i][c] -= dissolved_qty;
        conserve_qty[i][c] += dissolved_qty;
      }

    }
  }
}


/* ******************************************************************
* Computes source and sink terms and adds them to vector tcc.
* Returns mass rate for the tracer.
* The routine treats two cases of tcc with one and all components.
****************************************************************** */
void
Transport_ATS::AddSourceTerms_(double t0,
        double t1,
        Epetra_MultiVector& conserve_qty,
        int n0,
        int n1)
{
  int nsrcs = srcs_.size();

  for (int m = 0; m < nsrcs; m++) {
    srcs_[m]->Compute(t0, t1);
    std::vector<int> tcc_index = srcs_[m]->tcc_index();

    for (auto it = srcs_[m]->begin(); it != srcs_[m]->end(); ++it) {
      int c = it->first;
      std::vector<double>& values = it->second;

      if (c < conserve_qty.MyLength()) {
        if (srcs_[m]->getType() == DomainFunction_kind::COUPLING && n0 == 0) {
          AMANZI_ASSERT(values.size() == num_aqueous_ + 1);
          conserve_qty[num_aqueous_][c] += values[num_aqueous_];
        }

        for (int k = 0; k < tcc_index.size(); ++k) {
          int i = tcc_index[k];
          if (i < n0 || i > n1) continue;

          int imap = i;
          double value = mesh_->getCellVolume(c) * values[k];
          conserve_qty[imap][c] += (t1 - t0) * value;
        }
      }
    }
  }
}


} // namespace Transport
} // namespace Amanzi
