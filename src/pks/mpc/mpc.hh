/*
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

#ifndef PKS_MPC_MPC_HH_
#define PKS_MPC_MPC_HH_

#include <vector>
#include "boost/algorithm/string.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"

#include "State.hh"
#include "TreeVector.hh"


#include "PK.hh"
#include "PK_Factory.hh"

namespace Amanzi {

template <class PK_t>
class MPC : virtual public PK {

public:

  MPC(Teuchos::ParameterList& pk_tree,
      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
      const Teuchos::RCP<State>& S,
      const Teuchos::RCP<TreeVector>& solution)
      : PK(pk_tree, global_list, S, solution),
        global_list_(global_list),
        pk_tree_(pk_tree),
        pks_list_(Teuchos::sublist(global_list, "PKs"))        
  {
    // name the PK
    name_ = Keys::cleanPListName(pk_tree.name());

    // get my parameter list
    plist_ = Teuchos::sublist(pks_list_, name_);

    // set the verbose object list if need be
    if (plist_->isSublist(name_ + " verbose object")) {
      plist_->set("verbose object", plist_->sublist(name_ + " verbose object"));
      std::cout << "Overwriting VO with name = " << name_ << std::endl;
    } else {
      std::cout << "Default VO with name = " << name_ << std::endl;
    }


    // verbose object
    vo_ = Teuchos::rcp(new VerboseObject(solution->Comm(), name_, *plist_));
  }

  // Virtual destructor
  virtual ~MPC() {}

  // PK methods
  // -- setup
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- calls all sub-PK initialize() methods
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // transfer operators
  virtual void State_to_Solution(const Teuchos::RCP<State>& S,
          TreeVector& soln);
  virtual void Solution_to_State(TreeVector& soln,
          const Teuchos::RCP<State>& S);
  virtual void Solution_to_State(const TreeVector& soln,
          const Teuchos::RCP<State>& S);

  // -- loops over sub-PKs
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S);
  virtual bool ValidStep();
  virtual void ChangedSolutionPK(const Teuchos::Ptr<State>& S);
  
  // set States
  virtual void set_states(const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);


 protected:
  // this does not work.  fail etc
  void init_(const Teuchos::RCP<State>& S);

 protected:
  
  typedef std::vector<Teuchos::RCP<PK_t> > SubPKList;
  Teuchos::RCP<Teuchos::ParameterList> global_list_;
  Teuchos::ParameterList pk_tree_;
  Teuchos::RCP<Teuchos::ParameterList> pks_list_;

  SubPKList sub_pks_;

};


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::Setup(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->Setup(S);
  }
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::Initialize(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->Initialize(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their state_to_solution method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::State_to_Solution(const Teuchos::RCP<State>& S,
                            TreeVector& soln) {
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
    Teuchos::RCP<TreeVector> pk_soln = soln.SubVector(i);
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    sub_pks_[i]->State_to_Solution(S, *pk_soln);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their solution_to_state method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::Solution_to_State(TreeVector& soln,
                            const Teuchos::RCP<State>& S) {
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
    Teuchos::RCP<TreeVector> pk_soln = soln.SubVector(i);
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    sub_pks_[i]->Solution_to_State(*pk_soln, S);
  }
};

template <class PK_t>
void MPC<PK_t>::Solution_to_State(const TreeVector& soln,
                            const Teuchos::RCP<State>& S) {
  TreeVector* soln_nc_ptr = const_cast<TreeVector*>(&soln);
  Solution_to_State(*soln_nc_ptr, S);
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit_state method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "commiting step" << std::endl;
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->CommitStep(t_old, t_new, S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their calculate_diagnostics method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::CalculateDiagnostics(const Teuchos::RCP<State>& S) {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "calculating diagnostics" << std::endl;
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->CalculateDiagnostics(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their ValidStep() method
// -----------------------------------------------------------------------------
template <class PK_t>
bool MPC<PK_t>::ValidStep() {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Validating time step." << std::endl;

  bool valid(true);
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    bool is_valid = (*pk)->ValidStep();
    if (!is_valid) {
      valid = is_valid;
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Invalid time step, sub_pk: " << (*pk)->name()
                   << " is invalid." << std::endl;
    }
  }
  return valid;
};


// -----------------------------------------------------------------------------
// Marks sub-PKs as changed.
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::ChangedSolutionPK(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->ChangedSolutionPK(S.ptr());
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their set_states() methods
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::set_states(const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<State>& S_inter,
                     const Teuchos::RCP<State>& S_next) {
  //  PKDefaultBase::set_states(S, S_inter, S_next);
  S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;

  //  do the loop
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_states(S, S_inter, S_next);
  }
};

// protected constructor of subpks
template <class PK_t>
void MPC<PK_t>::init_(const Teuchos::RCP<State>& S)
{
  PKFactory pk_factory;
  Teuchos::Array<std::string> pk_order = plist_->get< Teuchos::Array<std::string> >("PKs order");

  int npks = pk_order.size();
  for (int i=0; i!=npks; ++i) {
    // create the solution vector
    Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(solution_->Comm()));
    solution_->PushBack(pk_soln);

    // create the PK
    std::string name_i = pk_order[i];
    Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(name_i, pk_tree_, global_list_, S, pk_soln);
    Teuchos::RCP<PK_t> pk = Teuchos::rcp_dynamic_cast<PK_t>(pk_notype, true); 
    sub_pks_.push_back(pk);
  }
};

} // close namespace Amanzi

#endif
