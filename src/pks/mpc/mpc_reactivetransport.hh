/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/
/*!

This MPC couples transport and chemistry through an operator split strategy.

Note that transport's primary variable is `"mole_fraction`", which is in units
of [mol C mol liquid^-1].  Chemistry is in "total_component_concentration", which
is in units of [mol-C L^-1].  Therefore, between steps, we convert between the
two.

`"PK type`" = `"reactive transport`"

.. _pk-reactive-transport-spec:
.. admonition:: pk-reactive-transport-spec

   * `"domain name`" ``[string]`` Domain of simulation

   KEYS:

   - `"molar density liquid`"


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

  virtual void parseParameterList() override;

  virtual void Setup() override;
  virtual void Initialize() override;
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

 protected:
  virtual void cast_sub_pks_();

 protected:
  Key domain_;
  Key tcc_key_;
  Key mol_frac_key_;
  Key mol_dens_key_;

  Teuchos::RCP<Teuchos::Time> alquimia_timer_;

  // storage for the component concentration intermediate values
  Teuchos::RCP<Transport::Transport_ATS> transport_pk_;

#ifdef ENABLE_ALQUIMIA
  Teuchos::RCP<AmanziChemistry::Alquimia_PK> chemistry_pk_;
#else
  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chemistry_pk_;
#endif

  // factory registration
  static RegisteredPKFactory<MPCReactiveTransport> reg_;
};

} // namespace Amanzi
