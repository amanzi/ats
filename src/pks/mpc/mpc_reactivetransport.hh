/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  Process kernel for coupling of Transport_PK and Chemistry_PK.
*/


#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "PK.hh"
#include "transport_ats.hh"
#include "Chemistry_PK.hh"
#include "weak_mpc.hh"

namespace Amanzi {

class MPCReactiveTransport : public WeakMPC {
 public:
  MPCReactiveTransport(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual double get_dt() override;

  virtual void parseParameterList() override;

  virtual void Setup() override;
  virtual void Initialize() override;
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

 protected:
  virtual void cast_sub_pks_();

 protected:
  bool chem_step_succeeded_;
  bool trans_step_succeeded_;

  Key domain_;
  Key tcc_key_;
  Key mol_dens_key_;

  Teuchos::RCP<Teuchos::Time> alquimia_timer_;

  // storage for the component concentration intermediate values
  Teuchos::RCP<Transport::Transport_ATS> transport_pk_;
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> chemistry_pk_;

  // factory registration
  static RegisteredPKFactory<MPCReactiveTransport> reg_;
};

} // namespace Amanzi
