/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

/*
  This is the mpc_pk component of the Amanzi code.

  PK for coupling of surface and subsurface transport PKs
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
  void modifyParameterList() override;
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
