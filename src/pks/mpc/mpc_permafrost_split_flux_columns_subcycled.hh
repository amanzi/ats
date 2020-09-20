/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

An operator-split permafrost coupler, based on flux.

solve:

(dTheta_s / dt)^* = div k_s grad (z+h)

then solve:

dTheta_s / dt = (dTheta_s / dt)^* + Q_ext + q_ss
dTheta / dt = div k (grad p + rho*g*\hat{z})
k  (grad p + rho*g*\hat{z}) |_s = q_ss

This effectively does an operator splitting on the surface flow equation, but
instead of the typical strateegy of passing pressure, passes the divergence of
lateral fluxes as a fixed source term.

This is the permafrost analog, so deals with energy as well in a similar
strategy.  In this case advection and diffusion of energy are handled in the
first solve:

(dE_s / dt)^* = div (  kappa_s grad T + hq )

then:

dE_s / dt = (dE_s / dt)^* + QE_ext + h * Q_ext + qE_ss + h * q_ss
dE / dt = div (  kappa grad T) + hq )
kappa grad T |_s = qE_ss


------------------------------------------------------------------------- */

#ifndef PKS_MPC_PERMAFROST_SPLIT_FLUX_COLUMNS_SUBCYCLED_HH_
#define PKS_MPC_PERMAFROST_SPLIT_FLUX_COLUMNS_SUBCYCLED_HH_

#include "PK.hh"
#include "mpc.hh"
#include "primary_variable_field_evaluator.hh"
#include "mpc_permafrost_split_flux_columns.hh"

namespace Amanzi {

class MPCPermafrostSplitFluxColumnsSubcycled : public MPCPermafrostSplitFluxColumns {

 public:

  MPCPermafrostSplitFluxColumnsSubcycled(Teuchos::ParameterList& FElist,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~MPCPermafrostSplitFluxColumnsSubcycled() = default;

  virtual double get_dt();
  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  virtual bool ValidStep();
  
  virtual void CommitStep(double t_old, double t_new,
                          const Teuchos::RCP<State>& S);
  
 private:
  std::string subcycled_timestep_type_;
  double  subcycled_timestep_target_, subcycled_target_time_;
  bool surface_star_subcycling_;
  // factory registration
  static RegisteredPKFactory<MPCPermafrostSplitFluxColumnsSubcycled> reg_;
  

};
} // close namespace Amanzi

#endif
