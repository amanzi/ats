/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Interface for EWC, a helper class that does projections and preconditioners in
energy/water-content space instead of temperature/pressure space, in the
subsurface.
------------------------------------------------------------------------- */

#ifndef MPC_DELEGATE_EWC_SURFACE_HH_
#define MPC_DELEGATE_EWC_SURFACE_HH_

#include "mpc_delegate_ewc.hh"

namespace Amanzi {

class MPCDelegateEWCSurface : public MPCDelegateEWC {
 public:
  MPCDelegateEWCSurface(Teuchos::ParameterList& plist, const Teuchos::RCP<State>& S);

 protected:
  virtual bool modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> up);
  virtual void precon_ewc_(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

 protected:
  double T_cutoff_;
};

} // namespace Amanzi


#endif
