/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

Globally implicit coupling solves all sub-PKs as a single system of equations.
This can be completely automated when all PKs are also :ref:`PK: BDF` PKs,
using a block-diagonal preconditioner where each diagonal block is provided by
its own sub-PK.

`"PK type`" = `"strong MPC`"

.. _pk-strong-mpc-spec:
.. admonition:: pk-strong-mpc-spec

   INCLUDES:

   - ``[mpc-spec]`` *Is a* :ref:`MPC`.
   - ``[pk-bdf-default-spec]`` *Is a* :ref:`PK: BDF`.

*/

#pragma once

#include <vector>

#include "mpc.hh"
#include "pk_bdf_default.hh"

namespace Amanzi {

// note this looks odd, but StrongMPC is both a MPC within a hierarchy of BDF
// PKs, but it also IS a BDF PK itself, in that it implements the BDF
// interface and can be implicitly integrated.
template<class PK_t>
class StrongMPC
  : public MPC<PK_t>
  , public ATS_Physics::PK_BDF_Default {
 public:
  StrongMPC(Teuchos::ParameterList& pk_list,
            const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
            const Teuchos::RCP<State>& S,
            const Teuchos::RCP<TreeVector>& solution);

  virtual void parseParameterList() override;

  //
  // These methods override methods in MPC<PK_t>
  // -----------------------------------------------------------------------------
  virtual void Setup() override;
  virtual void Initialize() override;

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void FailStep(double t_old, double t_new, const Tag& tag) override;

  //
  // These methods override methods in PK_BDF_Default
  // -----------------------------------------------------------------------------

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk FunctionalResidual().
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g) override;

  // -- enorm for the coupled system
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) override;

  // StrongMPC's preconditioner is, by default, just the block-diagonal
  // operator formed by placing the sub PK's preconditioners on the diagonal.
  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> Pu) override;

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.

  virtual void ChangedSolution(const Tag& S) override;
  virtual void ChangedSolution() override;

  // -- Admissibility of the solution.
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> u) override;

  // Is the step valid?
  virtual bool IsValid(const Teuchos::RCP<const TreeVector>& u) override;

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h,
                               Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) override;

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult ModifyCorrection(
    double h,
    Teuchos::RCP<const TreeVector> res,
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> du) override;

 protected:
  using MPC<PK_t>::sub_pks_;
  using MPC<PK_t>::global_list_;
  using MPC<PK_t>::pk_tree_;
  using MPC<PK_t>::pks_list_;

 private:
  // factory registration
  static RegisteredPKFactory<StrongMPC> reg_;
};

//
// Class implementation
//

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template<class PK_t>
StrongMPC<PK_t>::StrongMPC(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln)
  : PK(pk_tree, global_list, S, soln),
    MPC<PK_t>(pk_tree, global_list, S, soln),
    ATS_Physics::PK_BDF_Default(pk_tree, global_list, S, soln)
{
  MPC<PK_t>::init_(soln->Comm());
}


template<class PK_t>
void
StrongMPC<PK_t>::parseParameterList()
{
  // push on a parameter to indicate that sub-pks need not assemble their
  // operators, as we will do that here (or above here)
  auto pk_order = plist_->template get<Teuchos::Array<std::string>>("PKs order");
  for (const auto& pk_name : pk_order) {
    pks_list_->sublist(pk_name).set("strongly coupled PK", true);
  }

  MPC<PK_t>::parseParameterList();
}


// -----------------------------------------------------------------------------
// Setup
// -----------------------------------------------------------------------------
template<class PK_t>
void
StrongMPC<PK_t>::Setup()
{
  // push on a parameter to indicate that sub-pks need not assemble their
  // operators, as we will do that here (or above here)
  auto pk_order = plist_->get<Teuchos::Array<std::string>>("PKs order");
  for (const auto& pk_name : pk_order) {
    pks_list_->sublist(pk_name).set("strongly coupled PK", true);
  }

  MPC<PK_t>::Setup();
  ATS_Physics::PK_BDF_Default::Setup();
};


// -----------------------------------------------------------------------------
// Required unique initialize(), just calls both of its base class
// initialize() methods.
// -----------------------------------------------------------------------------
template<class PK_t>
void
StrongMPC<PK_t>::Initialize()
{
  // Just calls both subclass's initialize.  NOTE - order is important here --
  // MPC<PK_t> grabs the primary variables from each sub-PK and stuffs
  // them into the solution, which must be done prior to BDFBase initializing
  // the timestepper.

  // Initialize all sub PKs.
  MPC<PK_t>::Initialize();

  // Initialize my timestepper.
  ATS_Physics::PK_BDF_Default::Initialize();
};


// -----------------------------------------------------------------------------
// Calls both parts
// -----------------------------------------------------------------------------
template<class PK_t>
void
StrongMPC<PK_t>::CommitStep(double t_old, double t_new, const Tag& tag)
{
  MPC<PK_t>::CommitStep(t_old, t_new, tag);
  ATS_Physics::PK_BDF_Default::CommitStep(t_old, t_new, tag);
}

template<class PK_t>
void
StrongMPC<PK_t>::FailStep(double t_old, double t_new, const Tag& tag)
{
  MPC<PK_t>::FailStep(t_old, t_new, tag);
  //PK_BDF_Default::FailStep(t_old, t_new, tag);
}


// -----------------------------------------------------------------------------
// Compute the non-linear functional g = g(t,u,udot).
// -----------------------------------------------------------------------------
template<class PK_t>
void
StrongMPC<PK_t>::FunctionalResidual(double t_old,
                                    double t_new,
                                    Teuchos::RCP<const TreeVector> u_old,
                                    Teuchos::RCP<TreeVector> u_new,
                                    Teuchos::RCP<TreeVector> g)
{
  Solution_to_State(*u_new, tag_next_);

  // loop over sub-PKs
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the old solution sub-vector
    Teuchos::RCP<const TreeVector> pk_u_old(Teuchos::null);
    if (u_old != Teuchos::null) {
      pk_u_old = u_old->SubVector(i);
      if (pk_u_old == Teuchos::null) {
        Errors::Message message("MPC: vector structure does not match PK structure");
        Exceptions::amanzi_throw(message);
      }
    }

    // pull out the new solution sub-vector
    Teuchos::RCP<TreeVector> pk_u_new = u_new->SubVector(i);
    if (pk_u_new == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the residual sub-vector
    Teuchos::RCP<TreeVector> pk_g = g->SubVector(i);
    if (pk_g == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // fill the nonlinear function with each sub-PKs contribution
    sub_pks_[i]->FunctionalResidual(t_old, t_new, pk_u_old, pk_u_new, pk_g);
  }
};


// -----------------------------------------------------------------------------
// Applies preconditioner to u and returns the result in Pu.
// -----------------------------------------------------------------------------
template<class PK_t>
int
StrongMPC<PK_t>::ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu)
{
  // loop over sub-PKs
  int ierr = 0;
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the preconditioned u sub-vector
    Teuchos::RCP<TreeVector> pk_Pu = Pu->SubVector(i);
    if (pk_Pu == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // Fill the preconditioned u as the block-diagonal product using each sub-PK.
    int icur_err = sub_pks_[i]->ApplyPreconditioner(pk_u, pk_Pu);
    ierr += icur_err;
  }
  return ierr;
};


// -----------------------------------------------------------------------------
// Compute a norm on u-du and returns the result.
// For a Strong MPC, the enorm is just the max of the sub PKs enorms.
// -----------------------------------------------------------------------------
template<class PK_t>
double
StrongMPC<PK_t>::ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du)
{
  double norm = 0.0;

  // loop over sub-PKs
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // pull out the du sub-vector
    Teuchos::RCP<const TreeVector> pk_du = du->SubVector(i);
    if (pk_du == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // norm is the max of the sub-PK norms
    double tmp_norm = sub_pks_[i]->ErrorNorm(pk_u, pk_du);
    norm = std::max(norm, tmp_norm);
  }
  return norm;
};


// -----------------------------------------------------------------------------
// Update the preconditioner.
// -----------------------------------------------------------------------------
template<class PK_t>
void
StrongMPC<PK_t>::UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h)
{
  Solution_to_State(*up, tag_next_);

  // loop over sub-PKs
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the up sub-vector
    Teuchos::RCP<const TreeVector> pk_up = up->SubVector(i);
    if (pk_up == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    // update precons of each of the sub-PKs
    sub_pks_[i]->UpdatePreconditioner(t, pk_up, h);
  };
};


// -----------------------------------------------------------------------------
// Experimental approach -- calling this indicates that the time integration
// scheme is changing the value of the solution in state.
// -----------------------------------------------------------------------------
template<class PK_t>
void
StrongMPC<PK_t>::ChangedSolution(const Tag& tag)
{
  // loop over sub-PKs
  for (auto& pk : sub_pks_) pk->ChangedSolution(tag);
};


// -----------------------------------------------------------------------------
// Calling this indicates that the time integration scheme is changing
// the value of the solution in state.
// -----------------------------------------------------------------------------
template<class PK_t>
void
StrongMPC<PK_t>::ChangedSolution()
{
  // loop over sub-PKs
  for (auto& pk : sub_pks_) pk->ChangedSolution();
};

// -----------------------------------------------------------------------------
// Check admissibility of each sub-pk
// -----------------------------------------------------------------------------
template<class PK_t>
bool
StrongMPC<PK_t>::IsAdmissible(Teuchos::RCP<const TreeVector> u)
{
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    if (!sub_pks_[i]->IsAdmissible(pk_u)) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "PK " << sub_pks_[i]->name() << " is not admissible." << std::endl;
      return false;
    }
  }
  return true;
};


// -----------------------------------------------------------------------------
// Check validity of each sub-pk
// -----------------------------------------------------------------------------
template<class PK_t>
bool
StrongMPC<PK_t>::IsValid(const Teuchos::RCP<const TreeVector>& u)
{
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    if (!sub_pks_[i]->IsValid(pk_u)) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "PK " << sub_pks_[i]->name() << " is not admissible." << std::endl;
      return false;
    }
  }
  return true;
};


// -----------------------------------------------------------------------------
// Modify predictor from each sub pk.
// -----------------------------------------------------------------------------
template<class PK_t>
bool
StrongMPC<PK_t>::ModifyPredictor(double h,
                                 Teuchos::RCP<const TreeVector> u0,
                                 Teuchos::RCP<TreeVector> u)
{
  // loop over sub-PKs
  bool modified = false;
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u0 = u0->SubVector(i);
    Teuchos::RCP<TreeVector> pk_u = u->SubVector(i);
    if (pk_u == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    modified |= sub_pks_[i]->ModifyPredictor(h, pk_u0, pk_u);
  }
  return modified;
};


// -----------------------------------------------------------------------------
// Modify correction from each sub pk.
// -----------------------------------------------------------------------------
template<class PK_t>
AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
StrongMPC<PK_t>::ModifyCorrection(double h,
                                  Teuchos::RCP<const TreeVector> res,
                                  Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> du)
{
  // loop over sub-PKs
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult modified =
    AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    // pull out the u sub-vector
    Teuchos::RCP<const TreeVector> pk_u = u->SubVector(i);
    Teuchos::RCP<const TreeVector> pk_res = res->SubVector(i);
    Teuchos::RCP<TreeVector> pk_du = du->SubVector(i);

    if (pk_u == Teuchos::null || pk_du == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    modified = std::max(modified, sub_pks_[i]->ModifyCorrection(h, pk_res, pk_u, pk_du));
  }
  return modified;
};


} // namespace Amanzi
