/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/
/*!

This MPC couples surface and subsurface transport.  It is implemented as
described in `Molins et al WRR 2022 <https://doi.org/10.1029/2022WR032074>`_,
and deals with the conservation of advected fluxes as they are transported
laterally, infiltrate, etc, including wetting and drying of surface cells.

`"PK type`" = `"surface subsurface transport`"

.. _pk-surface-subsurface-transport-spec:
.. admonition:: pk-surface-subsurface-transport-spec

   * `"PKs order`" ``[Array(string)]`` Order must be {subsurface transport,
     surface transport}.

*/

#pragma once

#include "Teuchos_RCP.hpp"

//#include "pk_mpcsubcycled_ats.hh"
#include "weak_mpc.hh"
#include "PK.hh"
#include "transport_ats.hh"

namespace Amanzi {

class MPCCoupledTransport : public WeakMPC {
 public:
  MPCCoupledTransport(Teuchos::ParameterList& pk_tree,
                      const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  void parseParameterList() override;
  void Setup() override;
  int get_num_aqueous_component();

  // bug, see amanzi/ats#125
  // Once that is fixed, remove this advance in favor of the default implementation.
  bool AdvanceStep(double t_old, double t_new, bool reinit) override;


 protected:
  void SetupCouplingConditions_();

 protected:
  Teuchos::RCP<Transport::Transport_ATS> pk_ss_, pk_surf_;
  Key name_ss_, name_surf_;

  // factory registration
  static RegisteredPKFactory<MPCCoupledTransport> reg_;
};

} // namespace Amanzi
