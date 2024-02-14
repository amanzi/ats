/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Multi process coupler base class.
/*!

A multi process coupler is a PK (process kernel) which coordinates and couples
several PKs.  Each of these coordinated PKs may be MPCs themselves, or physical
PKs.  Note this does NOT provide a full implementation of PK -- it does not
supply the AdvanceStep() method.  Therefore this class cannot be instantiated, but
must be inherited by derived classes which finish supplying the functionality.
Instead, this provides the data structures and methods (which may be overridden
by derived classes) for managing multiple PKs.

Most of these methods simply loop through the coordinated PKs, calling their
respective methods.

.. _mpc-spec:
.. admonition:: mpc-spec

    * `"PKs order`" ``[Array(string)]`` Provide a specific order to the
      sub-PKs; most methods loop over all sub-PKs, and will call the sub-PK
      method in this order.

    INCLUDES:

    - ``[pk-spec]`` *Is a* PK_.

*/

#pragma once

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TreeVector.hh"

#include "PK.hh"
#include "PKFactory.hh"

namespace Amanzi {

template <class PK_t>
class MPC : virtual public PK {
 public:
  MPC(const Comm_ptr_type& comm,
      Teuchos::ParameterList& pk_tree,
      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
      const Teuchos::RCP<State>& S)
    : PK(comm, pk_tree, global_list, S),
      pk_tree_(pk_tree),
      global_list_(global_list),
      pks_list_(Teuchos::sublist(global_list, "PKs"))
  {}

  // PK methods
  // -- setup
  virtual void Setup() override;

  // -- calls all sub-PK initialize() methods
  virtual void Initialize() override;

  // -- setters/getters
  virtual void setTags(const Tag& current, const Tag& next) override;

  // additional getter
  virtual Teuchos::RCP<PK_t> get_subpk(int i);

  // -- PK timestepping
  // This is called after ALL PKs have successfully advanced their
  // steps, so information needed to back up can be overwritten.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // This is called if ANY PK has failed; do what is needed to back up for a
  // new attempt at the step.
  virtual void FailStep(double t_old, double t_new, const Tag& tag) override;

  // Calculate any diagnostics at S->time(), currently for visualization.
  virtual void CalculateDiagnostics(const Tag& tag) override;

  // Is the step valid?
  virtual bool ValidStep() override;

  // -- transfer operators
  virtual Teuchos::RCP<TreeVectorSpace> getSolutionSpace() const override;

  virtual void State_to_Solution(const Tag& tag, TreeVector& soln) override;

  // why are there two here?  and is a non-const one necessary?
  //  virtual void Solution_to_State(TreeVector& soln, const Teuchos::RCP<State>& S) = 0;
  virtual void Solution_to_State(const TreeVector& soln, const Tag& tag) override;

  // Tag the primary variable as changed in the DAG
  virtual void ChangedSolutionPK(const Tag& tag) override;

 protected:
  // constructs sub-pks
  void createSubPKs_(Comm_ptr_type comm = Teuchos::null);

 protected:
  Teuchos::ParameterList pk_tree_;
  Teuchos::RCP<Teuchos::ParameterList> global_list_;

  Teuchos::RCP<Teuchos::ParameterList> pks_list_;
  using SubPKList = std::vector<Teuchos::RCP<PK_t>>;
  SubPKList sub_pks_;
};


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_t>
void
MPC<PK_t>::Setup()
{
  for (auto& pk : sub_pks_) pk->Setup();
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class PK_t>
void
MPC<PK_t>::Initialize()
{
  for (auto& pk : sub_pks_) pk->Initialize();
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling setTags
// -----------------------------------------------------------------------------
template <class PK_t>
void
MPC<PK_t>::setTags(const Tag& current, const Tag& next)
{
  PK::setTags(current, next);
  for (auto& pk : sub_pks_) pk->setTags(current, next);
}


template <class PK_t>
Teuchos::RCP<TreeVectorSpace>
MPC<PK_t>::getSolutionSpace() const
{
  auto tvs = Teuchos::rcp(new TreeVectorSpace(comm_));
  for (const auto& sub_pk : sub_pks_) { tvs->PushBack(sub_pk->getSolutionSpace()); }
  return tvs;
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their state_to_solution method
// -----------------------------------------------------------------------------
template <class PK_t>
void
MPC<PK_t>::State_to_Solution(const Tag& tag, TreeVector& soln)
{
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    Teuchos::RCP<TreeVector> pk_soln = soln.getSubVector(i);
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    sub_pks_[i]->State_to_Solution(tag, *pk_soln);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their solution_to_state method
// -----------------------------------------------------------------------------
template <class PK_t>
void
MPC<PK_t>::Solution_to_State(const TreeVector& soln, const Tag& tag)
{
  for (std::size_t i = 0; i != sub_pks_.size(); ++i) {
    auto pk_soln = soln.getSubVector(i);
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    sub_pks_[i]->Solution_to_State(*pk_soln, tag);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit_state method
// -----------------------------------------------------------------------------
template <class PK_t>
void
MPC<PK_t>::CommitStep(double t_old, double t_new, const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "commiting step @ " << tag << std::endl;
  for (auto& pk : sub_pks_) { pk->CommitStep(t_old, t_new, tag); }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their fail step method
// -----------------------------------------------------------------------------
template <class PK_t>
void
MPC<PK_t>::FailStep(double t_old, double t_new, const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "failing step @ " << tag << std::endl;
  for (auto& pk : sub_pks_) pk->FailStep(t_old, t_new, tag);
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their calculate_diagnostics method
// -----------------------------------------------------------------------------
template <class PK_t>
void
MPC<PK_t>::CalculateDiagnostics(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "calculating diagnostics @ " << tag << std::endl;
  for (auto& pk : sub_pks_) pk->CalculateDiagnostics(tag);
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their ValidStep
// -----------------------------------------------------------------------------
template <class PK_t>
bool
MPC<PK_t>::ValidStep()
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "Validating time step." << std::endl;

  for (auto& pk : sub_pks_) {
    bool valid = pk->ValidStep();
    if (!valid) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Invalid time step, sub_pk: " << pk->getName() << " is invalid." << std::endl;
      return valid;
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
// Marks sub-PKs as changed.
// -----------------------------------------------------------------------------
template <class PK_t>
void
MPC<PK_t>::ChangedSolutionPK(const Tag& tag)
{
  for (auto& pk : sub_pks_) pk->ChangedSolutionPK(tag);
};


template <class PK_t>
Teuchos::RCP<PK_t>
MPC<PK_t>::get_subpk(int i)
{
  if (i >= sub_pks_.size()) {
    return Teuchos::null;
  } else {
    return sub_pks_.at(i);
  }
}


// protected constructor of subpks
template <class PK_t>
void
MPC<PK_t>::createSubPKs_(Comm_ptr_type comm)
{
  PKFactory pk_factory;
  auto pk_order = plist_->get<Teuchos::Array<std::string>>("PKs order");
  if (comm == Teuchos::null) comm = comm_;

  int npks = pk_order.size();
  for (int i = 0; i != npks; ++i) {
    // create the PK
    std::string name_i = pk_order[i];
    Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(name_i, comm, pk_tree_, global_list_, S_);
    Teuchos::RCP<PK_t> pk = Teuchos::rcp_dynamic_cast<PK_t>(pk_notype, true);
    sub_pks_.push_back(pk);
  }
};


} // namespace Amanzi
