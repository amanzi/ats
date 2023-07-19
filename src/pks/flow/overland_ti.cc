/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Epetra_FECrsMatrix.h"
#include "EpetraExt_RowMatrixOut.h"
#include "boost/math/special_functions/fpclassify.hpp"

#include "overland.hh"
#include "Op.hh"

namespace Amanzi {
namespace Flow {

#define DEBUG_FLAG 1
#define DEBUG_ICE_FLAG 0
#define DEBUG_RES_FLAG 0


// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void
OverlandFlow::FunctionalResidual(double t_old,
                                 double t_new,
                                 Teuchos::RCP<TreeVector> u_old,
                                 Teuchos::RCP<TreeVector> u_new,
                                 Teuchos::RCP<TreeVector> g)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  niter_++;

  // bookkeeping
  double h = t_new - t_old;
  AMANZI_ASSERT(std::abs(S_->get_time(tag_inter_) - t_old) < 1.e-4 * h);
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t_new) < 1.e-4 * h);

  Teuchos::RCP<CompositeVector> u = u_new->Data();

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old << " t1 = " << t_new << " h = " << h
               << std::endl;
#endif

  // unnecessary here if not debeugging, but doesn't hurt either
  S_next_->GetEvaluator(Keys::getKey(domain_, "pres_elev"))->HasFieldChanged(S_next_.ptr(), name_);

#if DEBUG_FLAG
  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames;
  vnames.push_back("z");
  vnames.push_back("h_old");
  vnames.push_back("h_new");
  vnames.push_back("h+z");

  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.push_back(S_inter_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "elevation")).ptr());
  vecs.push_back(S_inter_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "ponded_depth")).ptr());
  vecs.push_back(S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "ponded_depth")).ptr());
  vecs.push_back(S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "pres_elev")).ptr());

  db_->WriteVectors(vnames, vecs, true);
#endif

  // pointer-copy temperature into state and update any auxilary data
  Solution_to_State(*u_new, S_next_);

  // update boundary conditions
  bc_head_->Compute(t_new);
  bc_flux_->Compute(t_new);
  bc_seepage_head_->Compute(t_new);
  UpdateBoundaryConditions_(S_next_.ptr());

  // diffusion term, treated implicitly
  ApplyDiffusion_(S_next_.ptr(), res.ptr());

#if DEBUG_FLAG
  vnames.resize(2);
  vecs.resize(2);
  vnames[0] = "k_s";
  vnames[1] = "uw k_s";
  vecs[0] = S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "overland_conductivity")).ptr();
  vecs[1] =
    S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "upwind_overland_conductivity")).ptr();
  db_->WriteVectors(vnames, vecs, true);
  db_->WriteVector("res (diff)", res.ptr(), true);
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("res (acc)", res.ptr(), true);
#endif

  // add rhs load value
  AddSourceTerms_(res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("res (src)", res.ptr(), true);
#endif
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
int
OverlandFlow::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon application:" << std::endl;

#if DEBUG_FLAG
  db_->WriteVector("h_res", u->Data().ptr(), true);
#endif

  // apply the preconditioner
  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());

#if DEBUG_FLAG
  db_->WriteVector("PC*h_res (h-coords)", Pu->Data().ptr(), true);
#endif

  return (ierr > 0) ? 0 : 1;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void
OverlandFlow::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4 * t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, S_next_);

  // calculating the operator is done in 3 steps:
  // 1. Diffusion components

  // 1.a: Pre-assembly updates.
  // -- update boundary condition markers, which set the BC type
  bc_head_->Compute(t);
  bc_flux_->Compute(t);
  bc_seepage_head_->Compute(t);
  UpdateBoundaryConditions_(S_next_.ptr());

  // -- update the rel perm according to the boundary info and upwinding
  // -- scheme of choice
  UpdatePermeabilityData_(S_next_.ptr());
  if (jacobian_) UpdatePermeabilityDerivativeData_(S_next_.ptr());

  Teuchos::RCP<const CompositeVector> cond =
    S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "upwind_overland_conductivity"));
  Teuchos::RCP<const CompositeVector> dcond = Teuchos::null;

  if (jacobian_) {
    if (preconditioner_->RangeMap().HasComponent("face")) {
      Key dkey = Keys::getDerivKey(Keys::getKey(domain_, "upwind_overland_conductivity"), key_);
      dcond = S_next_->GetPtr<CompositeVector>(dkey);
    } else {
      Key dkey = Keys::getDerivKey(Keys::getKey(domain_, "overland_conductivity"), key_);
      dcond = S_next_->GetPtr<CompositeVector>(dkey);
    }
  }

  // 1.b: Create all local matrices.
  preconditioner_->Init();
  preconditioner_diff_->SetScalarCoefficient(cond, dcond);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);
  if (jacobian_) {
    Teuchos::RCP<const CompositeVector> pres_elev = Teuchos::null;
    Teuchos::RCP<CompositeVector> flux = Teuchos::null;
    if (preconditioner_->RangeMap().HasComponent("face")) {
      flux = S_next_->GetPtrW<CompositeVector>("surface-water_flux", name_);
      preconditioner_diff_->UpdateFlux(pres_elev.ptr(), flux.ptr());
    } else {
      S_next_->GetEvaluator(Keys::getKey(domain_, "pres_elev"))
        ->HasFieldChanged(S_next_.ptr(), name_);
      pres_elev = S_next_->GetPtrW<CompositeVector>(Keys::getKey(domain_, "pres_elev"));
    }
    preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), pres_elev.ptr());
  }

  // 2. Accumulation shift
  //    The desire is to keep this matrix invertible for pressures less than
  //    atmospheric.  To do that, we keep the accumulation derivative
  //    non-zero, calculating dWC_bar / dh_bar, where bar indicates (p -
  //    p_atm), not max(p - p_atm,0).  Note that this operator is in h
  //    coordinates, not p coordinates, as the diffusion operator is applied
  //    to h.
  //

  // -- update the accumulation derivatives
  const auto& cv = *S_next_->GetPtr<CompositeVector>("surface-cell_volume");
  CompositeVector dwc_dh(cv, INIT_MODE_COPY);
  dwc_dh.Scale(1. / h);
  preconditioner_acc_->AddAccumulationTerm(dwc_dh, "cell");

  preconditioner_diff_->ApplyBCs(true, true, true);
};


// -----------------------------------------------------------------------------
// Default enorm that uses an abs and rel tolerance to monitor convergence.
// -----------------------------------------------------------------------------
double
OverlandFlow::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> res)
{
  const Epetra_MultiVector& pd =
    *S_next_->GetPtr<CompositeVector>(key_)->ViewComponent("cell", true);
  const Epetra_MultiVector& cv =
    *S_next_->GetPtr<CompositeVector>(cell_vol_key_)->ViewComponent("cell", true);

  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "ENorm (Infnorm) of: " << key_ << ": " << std::endl;

  Teuchos::RCP<const CompositeVector> dvec = res->Data();
  double h = S_->get_time(tag_next_) - S_->get_time(tag_inter_);

  double enorm_val = 0.0;
  for (CompositeVector::name_iterator comp = dvec->begin(); comp != dvec->end(); ++comp) {
    double enorm_comp = 0.0;
    int enorm_loc = -1;
    const Epetra_MultiVector& dvec_v = *dvec->ViewComponent(*comp, false);

    if (*comp == std::string("cell")) {
      // error done relative to extensive, conserved quantity
      int ncells = dvec->size(*comp, false);
      for (unsigned int c = 0; c != ncells; ++c) {
        double enorm_c =
          std::abs(h * dvec_v[0][c]) / (atol_ * cv[0][c] + rtol_ * std::abs(pd[0][c]) * cv[0][c]);

        if (enorm_c > enorm_comp) {
          enorm_comp = enorm_c;
          enorm_loc = c;
        }
      }

    } else if (*comp == std::string("face")) {
      // error in flux -- relative to cell's extensive conserved quantity
      int nfaces = dvec->size(*comp, false);
      bool scaled_constraint =
        plist_->sublist("diffusion").get<bool>("scaled constraint equation", true);
      double constraint_scaling_cutoff =
        plist_->sublist("diffusion").get<double>("constraint equation scaling cutoff", 1.0);
      const Epetra_MultiVector& kr_f =
        *S_next_
           ->GetPtrW<CompositeVector>(
             Keys::getDerivKey(Keys::getKey(domain_, "upwind_overland_conductivity"), key_))
           ->ViewComponent("face", false);

      for (unsigned int f = 0; f != nfaces; ++f) {
        AmanziMesh::Entity_ID_List cells;
        cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::OWNED);
        double cv_min =
          cells.size() == 1 ? cv[0][cells[0]] : std::min(cv[0][cells[0]], cv[0][cells[1]]);
        double conserved_min = cells.size() == 1 ? pd[0][cells[0]] * cv[0][cells[0]] :
                                                   std::min(pd[0][cells[0]] * cv[0][cells[0]],
                                                            pd[0][cells[1]] * cv[0][cells[1]]);

        double enorm_f = fluxtol_ * h * std::abs(dvec_v[0][f]) /
                         (atol_ * cv_min + rtol_ * std::abs(conserved_min));
        if (scaled_constraint && (kr_f[0][f] < constraint_scaling_cutoff)) enorm_f *= kr_f[0][f];

        if (enorm_f > enorm_comp) {
          enorm_comp = enorm_f;
          enorm_loc = f;
        }
      }

    } else {
      Errors::Message msg;
      msg << "Unused error component \"" << *comp << "\" in conserved quantity error norm.";
      Exceptions::amanzi_throw(msg);
    }

    // Write out Inf norms too.
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      double infnorm(0.);
      dvec_v.NormInf(&infnorm);

      ENorm_t err;
      ENorm_t l_err;
      l_err.value = enorm_comp;
      l_err.gid = dvec_v.Map().GID(enorm_loc);

      int ierr = MPI_Allreduce(&l_err, &err, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
      AMANZI_ASSERT(!ierr);
      *vo_->os() << "  ENorm (" << *comp << ") = " << err.value << "[" << err.gid << "] ("
                 << infnorm << ")" << std::endl;
    }

    enorm_val = std::max(enorm_val, enorm_comp);
  }

  double enorm_val_l = enorm_val;
  int ierr = MPI_Allreduce(&enorm_val_l, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  AMANZI_ASSERT(!ierr);
  return enorm_val;
};


} // namespace Flow
} // namespace Amanzi
