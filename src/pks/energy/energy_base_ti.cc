/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "boost/math/special_functions/fpclassify.hpp"
#include "boost/test/floating_point_comparison.hpp"

#include "Debugger.hh"
#include "BoundaryFunction.hh"
#include "Evaluator.hh"
#include "energy_base.hh"
#include "Op.hh"

namespace Amanzi {
namespace Energy {

#define DEBUG_FLAG 1
#define MORE_DEBUG_FLAG 0

// EnergyBase is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void
EnergyBase::FunctionalResidual(double t_old,
                               double t_new,
                               Teuchos::RCP<TreeVector> u_old,
                               Teuchos::RCP<TreeVector> u_new,
                               Teuchos::RCP<TreeVector> g)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // increment, get timestep
  niter_++;
  double h = t_new - t_old;
  AMANZI_ASSERT(std::abs(S_->get_time(tag_current_) - t_old) < 1.e-4 * h);
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t_new) < 1.e-4 * h);

  // pointer-copy temperature into states and update any auxilary data
  Solution_to_State(*u_new, tag_next_);
  Teuchos::RCP<CompositeVector> u = u_new->Data();

#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "----------------------------------------------------------------" << std::endl
               << "Residual calculation: t0 = " << t_old << " t1 = " << t_new << " h = " << h
               << std::endl;

  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames{ "T_old", "T_new" };
  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.emplace_back(S_->GetPtr<CompositeVector>(key_, tag_current_).ptr());
  vecs.emplace_back(u.ptr());


  db_->WriteVectors(vnames, vecs, true);

  // vnames[0] = "sl"; vnames[1] = "si";
  // vecs[0] = S_next_->GetFieldData("saturation_liquid").ptr();
  // vecs[1] = S_next_->GetFieldData("saturation_ice").ptr();
  // db_->WriteVectors(vnames, vecs, false);
#endif

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);
  db_->WriteBoundaryConditions(bc_markers(), bc_values());

  // zero out residual
  Teuchos::RCP<CompositeVector> res = g->Data();
  res->PutScalar(0.0);

  // diffusion term, implicit
  ApplyDiffusion_(tag_next_, res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("K", S_->GetPtr<CompositeVector>(conductivity_key_, tag_next_).ptr(), true);
  db_->WriteVector("res (diff)", res.ptr(), true);
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  vnames = { "e_old", "e_new" };
  vecs = { S_->GetPtr<CompositeVector>(conserved_key_, tag_current_).ptr(),
           S_->GetPtr<CompositeVector>(conserved_key_, tag_next_).ptr() };
  db_->WriteVectors(vnames, vecs, true);
  db_->WriteVector("res (acc)", res.ptr());
#endif

  // advection term
  if (is_advection_term_) {
    if (implicit_advection_) {
      AddAdvection_(tag_next_, res.ptr(), true);
    } else {
      AddAdvection_(tag_current_, res.ptr(), true);
    }
#if DEBUG_FLAG
    db_->WriteVector("res (adv)", res.ptr(), true);
#endif
  }

  // source terms
  AddSources_(tag_next_, res.ptr());
#if DEBUG_FLAG
  db_->WriteVector("res (src)", res.ptr());
#endif
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
int
EnergyBase::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
#if DEBUG_FLAG
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon application:" << std::endl;
  db_->WriteVector("T_res", u->Data().ptr(), true);
#endif

  // apply the preconditioner
  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());

#if DEBUG_FLAG
  db_->WriteVector("PC*T_res", Pu->Data().ptr(), true);
#endif

  return (ierr > 0) ? 0 : 1;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void
EnergyBase::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon update at t = " << t << std::endl;

  // update state with the solution up.
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4 * t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, tag_next_);

  // update boundary conditions
  ComputeBoundaryConditions_(tag_next_);
  UpdateBoundaryConditions_(tag_next_);

  // div K_e grad u
  // force mass matrices to change
  if (!deform_key_.empty() &&
      S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_ + " precon"))
    preconditioner_diff_->SetTensorCoefficient(Teuchos::null);

  UpdateConductivityData_(tag_next_);
  if (jacobian_) UpdateConductivityDerivativeData_(tag_next_);

  // jacobian term
  Teuchos::RCP<const CompositeVector> dKdT = Teuchos::null;
  if (jacobian_) {
    if (!duw_conductivity_key_.empty()) {
      dKdT = S_->GetPtr<CompositeVector>(duw_conductivity_key_, tag_next_);
    } else {
      dKdT = S_->GetDerivativePtr<CompositeVector>(conductivity_key_, tag_next_, key_, tag_next_);
    }
  }

  // -- primary term
  Teuchos::RCP<const CompositeVector> conductivity =
    S_->GetPtr<CompositeVector>(uw_conductivity_key_, tag_next_);
  Teuchos::RCP<const CompositeVector> temp = S_->GetPtr<CompositeVector>(key_, tag_next_);

  // create local matrices
  preconditioner_->Init();
  preconditioner_diff_->SetScalarCoefficient(conductivity, dKdT);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, temp.ptr());
  preconditioner_diff_->ApplyBCs(true, true, true);

  // -- local matrices, Jacobian term
  if (jacobian_) {
    Teuchos::RCP<CompositeVector> flux =
      S_->GetPtrW<CompositeVector>(energy_flux_key_, tag_next_, name_);
    preconditioner_diff_->UpdateFlux(up->Data().ptr(), flux.ptr());
    preconditioner_diff_->UpdateMatricesNewtonCorrection(flux.ptr(), up->Data().ptr());
  }

  // -- update the accumulation derivatives, de/dT
  S_->GetEvaluator(conserved_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
  const auto& de_dT =
    *S_->GetDerivativePtr<CompositeVector>(conserved_key_, tag_next_, key_, tag_next_)
       ->ViewComponent("cell", false);
  CompositeVector acc(S_->GetPtr<CompositeVector>(conserved_key_, tag_next_)->Map());
  auto& acc_c = *acc.ViewComponent("cell", false);


#if DEBUG_FLAG
  db_->WriteVector(
    "    de_dT",
    S_->GetDerivativePtr<CompositeVector>(conserved_key_, tag_next_, key_, tag_next_).ptr());
#endif

  unsigned int ncells = de_dT.MyLength();
  if (coupled_to_subsurface_via_temp_ || coupled_to_subsurface_via_flux_) {
    // do not add in de/dT if the height is 0
    const auto& pres = *S_->Get<CompositeVector>(Keys::getKey(domain_, "pressure"), tag_next_)
                          .ViewComponent("cell", false);
    const double& p_atm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    for (unsigned int c = 0; c != ncells; ++c) {
      acc_c[0][c] = pres[0][c] >= p_atm ? de_dT[0][c] / h : 0.;
    }
  } else {
    if (precon_used_) {
      for (unsigned int c = 0; c != ncells; ++c) {
        //      AMANZI_ASSERT(de_dT[0][c] > 1.e-10);
        // ?? Not using e_bar anymore apparently, though I didn't think we were ever.  Need a nonzero here to ensure not singlar.
        acc_c[0][c] = std::max(de_dT[0][c], 1.e-12) / h;
      }
    } else {
      if (decoupled_from_subsurface_) {
        const auto& uf_c =
          *S_->Get<CompositeVector>(uf_key_, tag_next_).ViewComponent("cell", false);
        for (unsigned int c = 0; c != ncells; ++c) {
          acc_c[0][c] = std::max(de_dT[0][c] / h, 1.e-1 * uf_c[0][c]) + 1.e-6;
        }

      } else {
        for (unsigned int c = 0; c != ncells; ++c) {
          //      AMANZI_ASSERT(de_dT[0][c] > 1.e-10);
          // ?? Not using e_bar anymore apparently, though I didn't think we were ever.  Need a nonzero here to ensure not singlar.
          // apply a diagonal shift manually for coupled problems
          acc_c[0][c] = de_dT[0][c] / h + 1.e-6;
        }
      }
    }
  }

  preconditioner_acc_->AddAccumulationTerm(acc, "cell");

  // -- update preconditioner with source term derivatives if needed
  AddSourcesToPrecon_(h);

  // update with advection terms
  if (is_advection_term_) {
    if (implicit_advection_ && implicit_advection_in_pc_) {
      Teuchos::RCP<const CompositeVector> water_flux =
        S_->GetPtr<CompositeVector>(flux_key_, tag_next_);
      S_->GetEvaluator(enthalpy_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
      Teuchos::RCP<const CompositeVector> dhdT =
        S_->GetDerivativePtr<CompositeVector>(enthalpy_key_, tag_next_, key_, tag_next_);
      preconditioner_adv_->Setup(*water_flux);
      preconditioner_adv_->SetBCs(bc_adv_, bc_adv_);
      preconditioner_adv_->UpdateMatrices(water_flux.ptr(), dhdT.ptr());
      preconditioner_adv_->ApplyBCs(false, true, false);
    }
  }

  // Apply boundary conditions.
  preconditioner_diff_->ApplyBCs(true, true, true);
};

// -----------------------------------------------------------------------------
// Default enorm that uses an abs and rel tolerance to monitor convergence.
// -----------------------------------------------------------------------------
double
EnergyBase::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> res)
{
  // Abs tol based on old conserved quantity -- we know these have been vetted
  // at some level whereas the new quantity is some iterate, and may be
  // anything from negative to overflow.

  //S_->GetEvaluator(conserved_key_, tag_current_).Update(*S_, name());
  // not used ?? jjb
  const Epetra_MultiVector& conserved =
    *S_->Get<CompositeVector>(conserved_key_, tag_current_).ViewComponent("cell", true);
  //S_->GetEvaluator(wc_key_, tag_current_).Update(*S_, name());
  const Epetra_MultiVector& wc =
    *S_->Get<CompositeVector>(wc_key_, tag_current_).ViewComponent("cell", true);
  const Epetra_MultiVector& cv =
    *S_->Get<CompositeVector>(cell_vol_key_, tag_next_).ViewComponent("cell", true);

  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_MEDIUM))
    *vo_->os() << "ENorm (Infnorm) of: " << conserved_key_ << ": " << std::endl;

  Teuchos::RCP<const CompositeVector> dvec = res->Data();
  double h = S_->get_time(tag_next_) - S_->get_time(tag_current_);

  Teuchos::RCP<const Comm_type> comm_p = mesh_->getComm();
  Teuchos::RCP<const MpiComm_type> mpi_comm_p =
    Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm_p);
  const MPI_Comm& comm = mpi_comm_p->Comm();

  double enorm_val = 0.0;
  for (CompositeVector::name_iterator comp = dvec->begin(); comp != dvec->end(); ++comp) {
    double enorm_comp = 0.0;
    int enorm_loc = -1;
    const Epetra_MultiVector& dvec_v = *dvec->ViewComponent(*comp, false);

    if (*comp == std::string("cell")) {
      // error done in two parts, relative to mass but absolute in
      // energy since it doesn't make much sense to be relative to
      // energy
      int ncells = dvec->size(*comp, false);
      for (unsigned int c = 0; c != ncells; ++c) {
        double mass = std::max(mass_atol_, wc[0][c] / cv[0][c]);
        double energy = mass * atol_ + soil_atol_;
        double enorm_c = std::abs(h * dvec_v[0][c]) / (energy * cv[0][c]);

        if (enorm_c > enorm_comp) {
          enorm_comp = enorm_c;
          enorm_loc = c;
        }
      }

    } else if (*comp == std::string("face")) {
      // error in flux -- relative to cell's extensive conserved quantity
      int nfaces = dvec->size(*comp, false);

      for (unsigned int f = 0; f != nfaces; ++f) {
        auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::OWNED);
        double cv_min =
          cells.size() == 1 ? cv[0][cells[0]] : std::min(cv[0][cells[0]], cv[0][cells[1]]);
        double mass_min = cells.size() == 1 ? wc[0][cells[0]] / cv[0][cells[0]] :
                                              std::min(wc[0][cells[0]] / cv[0][cells[0]],
                                                       wc[0][cells[1]] / cv[0][cells[1]]);
        mass_min = std::max(mass_min, mass_atol_);

        double energy = mass_min * atol_ + soil_atol_;
        double enorm_f = fluxtol_ * h * std::abs(dvec_v[0][f]) / (energy * cv_min);

        if (enorm_f > enorm_comp) {
          enorm_comp = enorm_f;
          enorm_loc = f;
        }
      }

    } else {
      // boundary face components had better be effectively identically 0
      double norm;
      dvec_v.Norm2(&norm);
      AMANZI_ASSERT(norm < 1.e-15);
    }


    // Write out Inf norms too.
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      double infnorm(0.);
      dvec_v.NormInf(&infnorm);

      ENorm_t err;
      ENorm_t l_err;
      l_err.value = enorm_comp;
      l_err.gid = dvec_v.Map().GID(enorm_loc);

      int ierr;
      ierr = MPI_Allreduce(&l_err, &err, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);
      AMANZI_ASSERT(!ierr);
      *vo_->os() << "  ENorm (" << *comp << ") = " << err.value << "[" << err.gid << "] ("
                 << infnorm << ")" << std::endl;
    }

    enorm_val = std::max(enorm_val, enorm_comp);
  }

  double enorm_val_l = enorm_val;

  int ierr;
  ierr = MPI_Allreduce(&enorm_val_l, &enorm_val, 1, MPI_DOUBLE, MPI_MAX, comm);
  AMANZI_ASSERT(!ierr);
  return enorm_val;
};


} // namespace Energy
} // namespace Amanzi
