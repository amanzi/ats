/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "boost/math/special_functions/fpclassify.hpp"
#include "Op.hh"
#include "snow_distribution.hh"

#define DEBUG_FLAG 1

namespace Amanzi {
namespace Flow {

// Overland is a BDFFnBase
// -----------------------------------------------------------------------------
// computes the non-linear functional g = g(t,u,udot)
// -----------------------------------------------------------------------------
void
SnowDistribution::FunctionalResidual(double t_old,
                                     double t_new,
                                     Teuchos::RCP<TreeVector> u_old,
                                     Teuchos::RCP<TreeVector> u_new,
                                     Teuchos::RCP<TreeVector> g)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();

  // bookkeeping
  double h = t_new - t_old;
  bool asserting = false;
  if (!(std::abs(S_->get_time(tag_current_) - t_old) < 1.e-4 * h)) {
    Errors::Message message;
    message << "SnowDistribution PK: ASSERT!: tag_current_ = " << S_->get_time(tag_current_)
            << ", t_old = " << t_old
            << " --> diff = " << std::abs(S_->get_time(tag_current_) - t_old) << " relative to "
            << h * 1.e-4
            << "  Maybe you checkpoint-restarted from a checkpoint file that was not on an even "
               "day?  This breaks the snow distribution!";
    Exceptions::amanzi_throw(message);
  }
  if (!(std::abs(S_->get_time(tag_next_) - t_new) < 1.e-4 * h)) {
    Errors::Message message;
    message << "SnowDistribution PK: ASSERT!: tag_next_ = " << S_->get_time(tag_next_)
            << ", t_new = " << t_new << " --> diff = " << std::abs(S_->get_time(tag_next_) - t_new)
            << " relative to " << h * 1.e-4
            << "  Maybe you checkpoint-restarted from a checkpoint file that was not on an even "
               "day?  This breaks the snow distribution!";
    Exceptions::amanzi_throw(message);
  }


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
  S_->GetEvaluator(potential_key_, tag_next_).Update(*S_, name_);

#if DEBUG_FLAG
  // dump u_old, u_new
  db_->WriteCellInfo(true);
  std::vector<std::string> vnames{ "precip_new", "potential" };
  std::vector<Teuchos::Ptr<const CompositeVector>> vecs;
  vecs.push_back(u.ptr());
  vecs.push_back(S_->GetPtr<CompositeVector>(potential_key_, tag_next_).ptr());

  db_->WriteVectors(vnames, vecs, true);
#endif

  // pointer-copy temperature into state and update any auxilary data
  Solution_to_State(*u_new, tag_next_);

  // diffusion term, treated implicitly
  ApplyDiffusion_(tag_next_, res.ptr());

#if DEBUG_FLAG
  db_->WriteVector("k_s", S_->GetPtrW<CompositeVector>(uw_cond_key_, tag_next_, name_).ptr(), true);
  db_->WriteVector("res (post diffusion)", res.ptr(), true);
#endif

  // accumulation term
  AddAccumulation_(res.ptr());
#if DEBUG_FLAG
  if (vo_->os_OK(Teuchos::VERB_HIGH)) {
    std::vector<double> time(1, S_->get_time(tag_next_));
    double precip = (*precip_func_)(time);
    *vo_->os() << " precip = " << precip << std::endl;
  }
  db_->WriteVector("res (post accumulation)", res.ptr(), true);
#endif
};


// -----------------------------------------------------------------------------
// Apply the preconditioner to u and return the result in Pu.
// -----------------------------------------------------------------------------
int
SnowDistribution::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH)) *vo_->os() << "Precon application:" << std::endl;

#if DEBUG_FLAG
  db_->WriteVector("h_res", u->Data().ptr(), true);
#endif

  // apply the preconditioner
  int ierr = preconditioner_->ApplyInverse(*u->Data(), *Pu->Data());
  Pu->Data()->Scale(1. / (10 * dt_factor_));


#if DEBUG_FLAG
  db_->WriteVector("PC*h_res", Pu->Data().ptr(), true);
#endif

  return (ierr > 0) ? 0 : 1;
};


// -----------------------------------------------------------------------------
// Update the preconditioner at time t and u = up
// -----------------------------------------------------------------------------
void
SnowDistribution::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  // VerboseObject stuff.
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Precon update at t = " << t << std::endl;


  // update state with the solution up.
  AMANZI_ASSERT(std::abs(S_->get_time(tag_next_) - t) <= 1.e-4 * t);
  PK_PhysicalBDF_Default::Solution_to_State(*up, tag_next_);

  // update the rel perm according to the scheme of choice
  UpdatePermeabilityData_(tag_next_);

  // primary variable
  Teuchos::RCP<const CompositeVector> cond = S_->GetPtr<CompositeVector>(uw_cond_key_, tag_next_);

  // Jacobian
  // playing it fast and loose.... --etc
  S_->GetEvaluator(cond_key_, tag_next_).UpdateDerivative(*S_, name_, key_, tag_next_);
  Teuchos::RCP<CompositeVector> dcond =
    S_->GetDerivativePtrW<CompositeVector>(cond_key_, tag_next_, key_, tag_next_, cond_key_);
  // NOTE: this scaling of dt is wrong, but keeps consistent with the diffusion derivatives
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);
  dcond->Scale(1. / (10 * dt_factor_));

  // calculating the operator is done in 3 steps:
  // 1. Create all local matrices.
  preconditioner_->Init();
  preconditioner_diff_->SetScalarCoefficient(cond, dcond);
  preconditioner_diff_->UpdateMatrices(Teuchos::null, Teuchos::null);

  S_->GetEvaluator(potential_key_, tag_next_).Update(*S_, name_);
  Teuchos::RCP<const CompositeVector> potential =
    S_->GetPtr<CompositeVector>(potential_key_, tag_next_);
  preconditioner_diff_->UpdateMatricesNewtonCorrection(Teuchos::null, potential.ptr());

  // 2.b Update local matrices diagonal with the accumulation terms.
  // -- update the accumulation derivatives

  const auto& cv = *S_->GetPtr<CompositeVector>(cv_key_, tag_next_);
  CompositeVector cv_times_10(cv);
  cv_times_10.Scale(10);
  preconditioner_acc_->AddAccumulationTerm(cv_times_10, "cell");

  preconditioner_diff_->ApplyBCs(true, true, true);
};

double
SnowDistribution::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  Teuchos::RCP<const CompositeVector> res = du->Data();
  const Epetra_MultiVector& res_c = *res->ViewComponent("cell", false);
  const Epetra_MultiVector& precip_c = *u->Data()->ViewComponent("cell", false);

  const Epetra_MultiVector& cv =
    *S_->Get<CompositeVector>(cv_key_, tag_next_).ViewComponent("cell", false);
  double dt = S_->get_time(tag_next_) - S_->get_time(tag_current_);
  std::vector<double> time(1, S_->get_time(tag_next_));

  // Cell error is based upon error in mass conservation
  double Qe = (*precip_func_)(time);
  double enorm_cell(0.);
  int bad_cell = -1;
  unsigned int ncells = res_c.MyLength();
  for (unsigned int c = 0; c != ncells; ++c) {
    double tmp =
      std::abs(res_c[0][c] * dt) / (atol_ * .1 * cv[0][c] + rtol_ * 10 * dt * precip_c[0][c]);
    if (tmp > enorm_cell) {
      enorm_cell = tmp;
      bad_cell = c;
    }
  }

  // Write out Inf norms too.
  if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
    double infnorm_c(0.);
    res_c.NormInf(&infnorm_c);

    ENorm_t err_c;
    ENorm_t l_err_c;
    l_err_c.value = enorm_cell;
    l_err_c.gid = res_c.Map().GID(bad_cell);

    MPI_Allreduce(&l_err_c, &err_c, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    *vo_->os() << "ENorm (cells) = " << err_c.value << "[" << err_c.gid << "] (" << infnorm_c << ")"
               << std::endl;
  }

  double buf = enorm_cell;
  MPI_Allreduce(&buf, &enorm_cell, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return enorm_cell;
};

bool
SnowDistribution::ModifyPredictor(double h,
                                  Teuchos::RCP<const TreeVector> u0,
                                  Teuchos::RCP<TreeVector> u)
{
  return false;
}


// Advance PK from time t_old to time t_new. True value of the last
// parameter indicates drastic change of boundary and/or source terms
// that may need PK's attention.
//
//  ALL SORTS OF FRAGILITY and UGLINESS HERE!
//  DO NOT USE THIS OUT IN THE WILD!
//
//  1. this MUST go first
//  2. it must be PERFECT NON_OVERLAPPING with everything else.  I'm not
//     sure exactly what that means.  Something like, nothing that this PK
//     writes can be read by anything else, except for the precip snow at
//     the end?  But it should be ok?
//  3. Exatrapolating in the timestepper should break things, so don't.
//  4. set: pk's distribution time, potential's dt factor
bool
SnowDistribution::AdvanceStep(double t_old, double t_new, bool reinit)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (t_new <= my_next_time_) {
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 << std::endl
                 << "BIG STEP still good!" << std::endl
                 << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 << std::endl;
    return false;
  }

  Teuchos::OSTab out = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_HIGH))
    *vo_->os() << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
               << std::endl
               << "BIG STEP Advancing: t0 = " << t_old << " t1 = " << t_old + dt_factor_
               << " h = " << dt_factor_ << std::endl
               << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
               << std::endl;

  std::vector<double> time(1, t_old);
  double Ps_old = (*precip_func_)(time);
  time[0] = t_old + dt_factor_;
  double Ps_new = (*precip_func_)(time);
  double Ps_mean = (Ps_new + Ps_old) / 2.;
  S_->GetW<CompositeVector>(key_, tag_current_, name_).PutScalar(Ps_mean);
  S_->GetW<CompositeVector>(key_, tag_next_, name_).PutScalar(Ps_mean);

  double my_dt = -1;
  double my_t_old = t_old;
  double my_t_new = t_old;

  Tag tag_subcycle_current(tag_current_);
  Tag tag_subcycle_next(tag_next_);
  S_->set_time(tag_subcycle_current, t_old);
  while (my_t_old < t_old + dt_factor_) {
    my_dt = PK_PhysicalBDF_Default::get_dt();
    my_t_new = std::min(my_t_old + my_dt, t_old + dt_factor_);
    S_->set_time(tag_subcycle_next, my_t_new);

    bool failed = PK_PhysicalBDF_Default::AdvanceStep(my_t_old, my_t_new, false);

    if (failed) {
      S_->set_time(tag_subcycle_next, S_->get_time(tag_subcycle_current));
    } else {
      S_->set_time(tag_subcycle_current, S_->get_time(tag_subcycle_next));
      PK_PhysicalBDF_Default::CommitStep(my_t_old, my_t_new, tag_subcycle_next);
      my_t_old = my_t_new;
    }
  }

  my_next_time_ = t_old + dt_factor_;

  // commit the precip to the OLD time as well -- this ensures that
  // even if a coupled PK fails at any point in the coming
  // distribution time, we keep the new value.
  *S_->GetPtrW<CompositeVector>(key_, tag_current_, name_) =
    *S_->GetPtr<CompositeVector>(key_, tag_next_);

  // clean up
  S_->set_time(tag_next_, t_new);
  S_->set_time(tag_current_, t_old);


  return false;
}


} // namespace Flow
} // namespace Amanzi
