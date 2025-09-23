/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*!

This couples integrated (surface + subsurface) transport MPC with an integrated
chemistry MPC.

Transport must be tightly coupled in order to correctly pass surface-subsurface
advective fluxes between the two domains.  Chemistry is just a weak MPC.  The
two are then coupled to form an integrated reactive-transport PK.

This performs operatoring splitting between transport and chemistry in each
domain.  Note that transport's primary variable is "molar_fraction", which is
in units of [mol-C mol-H2O^-1].  Chemistry is in
"total_component_concentration", which is in units of [mol-C L^-1].  Therefore,
between steps, we convert between the two.

`"PK type`" = `"surface subsurface reactive transport`"

.. _pk-surface-subsurface-reactive-transport-spec:
.. admonition:: pk-surface-subsurface-reactive-transport-spec

   * `"PKs order`" ``[Array(string)]`` Order must be {chemistry_mpc,
     transport_mpc}.  The chemistry MPC is likely just a :ref:`Weak MPC`
     coupling to chemistry PKs, while the transport MPC is likely a
     :ref:`Integrated Transport` MPC.

   KEYS:

   - `"molar density liquid`"
   - `"surface molar density liquid`"


*/


#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "PK.hh"
#include "mpc_coupled_transport.hh"
#include "transport_ats.hh"
#include "Chemistry_PK.hh"
#include "weak_mpc.hh"

namespace Amanzi {

class MPCCoupledReactiveTransport : public WeakMPC {
 public:
  MPCCoupledReactiveTransport(Teuchos::ParameterList& pk_tree,
                              const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                              const Teuchos::RCP<State>& S,
                              const Teuchos::RCP<TreeVector>& soln);


  void parseParameterList() override;

  // PK methods
  virtual void Setup() override;
  virtual void Initialize() override;
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

 protected:
  virtual void cast_sub_pks_();

 protected:
  Key domain_, domain_surf_;
  Key tcc_key_, tcc_surf_key_;
  Key mol_frac_key_, mol_frac_surf_key_;
  Key mol_dens_key_, mol_dens_surf_key_;

  // storage for the component concentration intermediate values
  Teuchos::RCP<MPCCoupledTransport> coupled_transport_pk_;
  Teuchos::RCP<WeakMPC> coupled_chemistry_pk_;

  Teuchos::RCP<ATS_Physics::Transport::Transport_ATS> transport_pk_;
  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chemistry_pk_;
  Teuchos::RCP<ATS_Physics::Transport::Transport_ATS> transport_pk_surf_;
  Teuchos::RCP<AmanziChemistry::Chemistry_PK> chemistry_pk_surf_;

  // factory registration
  static RegisteredPKFactory<MPCCoupledReactiveTransport> reg_;
};

} // namespace Amanzi
