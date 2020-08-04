/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! A DomainSet coupler, weakly couples PKs on domains of the same structure.
/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!
  
 */

#pragma once

#include "Key.hh"
#include "PK.hh"
#include "mpc.hh"

namespace Amanzi {

class DomainSetMPC : public MPC<PK> {
 public:
  DomainSetMPC(Teuchos::ParameterList& pk_tree,
               const Teuchos::RCP<Teuchos::ParameterList>& global_list,
               const Teuchos::RCP<State>& S,
               const Teuchos::RCP<TreeVector>& solution);

  virtual ~DomainSetMPC() = default;
  
  // PK methods
  virtual double get_dt();
  virtual void set_dt(double dt);
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

 protected:
  std::string pks_set_;

 private:
  // factory registration
  static RegisteredPKFactory<DomainSetMPC> reg_;
};

} // namspace
