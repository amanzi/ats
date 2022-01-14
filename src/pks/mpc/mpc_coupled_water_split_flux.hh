/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

An operator-split integrated hydrology coupler, based on flux.

solve:

(dTheta_s / dt)^* = div k_s grad (z+h)

then solve:

dTheta_s / dt = (dTheta_s / dt)^* + Q_ext + q_ss
dTheta / dt = div k (grad p + rho*g*\hat{z})
k  (grad p + rho*g*\hat{z}) |_s = q_ss

This effectively does an operator splitting on the surface flow equation, but
instead of the typical strateegy of passing pressure, passes the divergence of
lateral fluxes as a fixed source term.


------------------------------------------------------------------------- */

#ifndef PKS_MPC_COUPLED_WATER_SPLIT_FLUX_HH_
#define PKS_MPC_COUPLED_WATER_SPLIT_FLUX_HH_

#include "PK.hh"
#include "mpc.hh"
#include "EvaluatorPrimary.hh"

namespace Amanzi {

class MPCCoupledWaterSplitFlux : public MPC<PK> {

 public:

  MPCCoupledWaterSplitFlux(Teuchos::ParameterList& FElist,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~MPCCoupledWaterSplitFlux() = default;

  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();

  virtual void Initialize(const Teuchos::Ptr<State>& S);
  virtual void Setup(const Teuchos::Ptr<State>& S);
  
  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  virtual void set_dt(double dt);

  virtual void CommitStep(double t_old, double t_new,
                          const Teuchos::RCP<State>& S);
  
  virtual void CopyPrimaryToStar(const Teuchos::Ptr<const State>& S,
          const Teuchos::Ptr<State>& S_star);
  virtual void CopyStarToPrimary(double dt);

 protected:
  Key primary_variable_;
  Key primary_variable_star_;
  Key conserved_variable_star_;
  Key lateral_flow_source_;
  Key cv_key_;
  Teuchos::RCP<EvaluatorPrimary> eval_pvfe_;
  
 private:
  // factory registration
  static RegisteredPKFactory<MPCCoupledWaterSplitFlux> reg_;


};
} // close namespace Amanzi

#endif
