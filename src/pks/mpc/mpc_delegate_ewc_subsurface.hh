/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Globalization for nonlinearity associated with phase change and latent heat.

/*!

Interface for EWC, a helper class that does projections and preconditioners in
energy/water-content space instead of temperature/pressure space.

*/


#ifndef MPC_DELEGATE_EWC_SUBSURFACE_HH_
#define MPC_DELEGATE_EWC_SUBSURFACE_HH_

#include "mpc_delegate_ewc.hh"

namespace Amanzi {

class MPCDelegateEWCSubsurface : public MPCDelegateEWC {
 public:
  MPCDelegateEWCSubsurface(Teuchos::ParameterList& plist) :
      MPCDelegateEWC(plist) {}  

 protected:
  virtual bool modify_predictor_smart_ewc_(double h, Teuchos::RCP<TreeVector> up);
  virtual void precon_ewc_(Teuchos::RCP<const TreeVector> u,
                             Teuchos::RCP<TreeVector> Pu);

};

} // namespace


#endif
