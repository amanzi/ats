/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! An operator-split permafrost coupler, splitting overland flow from subsurface.

/*!
solve:

(dTheta_s / dt)^* = div k_s grad (z+h)

then solve:

dTheta_s / dt = (dTheta_s / dt)^* + Q_ext + q_ss
dTheta / dt = div k (grad p + rho*g*\hat{z})
k  (grad p + rho*g*\hat{z}) |_s = q_ss

This effectively does an operator splitting on the surface flow equation,
passing some combination of pressure and divergence of fluxes to the
subsurface.

This is the permafrost analog, so deals with energy as well in a similar
strategy.  In this case advection and diffusion of energy are handled in the
first solve:

(dE_s / dt)^* = div (  kappa_s grad T + hq )

then:

dE_s / dt = (dE_s / dt)^* + QE_ext + h * Q_ext + qE_ss + h * q_ss
dE / dt = div (  kappa grad T) + hq )
kappa grad T |_s = qE_ss

Note that this can be used with either a 3D subsurface solve, by setting the
2nd sub-PK to be a 3D permafrost MPC, or a bunch of columns, but setting the
2nd sub-PK to be a DomainSetMPC.


.. _mpc-permafrost-split-flux-spec
.. admonition:: mpc-permafrost-split-flux-spec

   * `"domain name`" ``[string]`` The subsurface domain, e.g. "domain" (for a
     3D subsurface ) or "column:*" (for the intermediate scale model.

   * `"star domain name`" ``[string]`` The surface domain, typically
     `"surface_star`" by convention.

   * `"coupling type`" ``[string]`` **hybrid** One of: `"pressure`" (pass the
     pressure field when coupling flow in the operator splitting), `"flux`"
     (pass the divergence of fluxes as a source), or `"hybrid`" a mixture of
     the two that seems the most robust.

   IF
   * `"subcycle subdomains`" ``[bool]`` **false** If true, subcycle surface_star
     system.  Typically this is paired with the `"subcycle subdomains`"
     option in the DomainSetMPC for columns in the ISM.

   THEN
   * `"subcycling target time step [s]`" ``[double]`` Step size to target or
     subcycling.  Note that this may be adjusted to meet events, etc.

   * `"minimum subcycled time step [s]`" ``[double]`` **1e-4** Errors out if
     any subdomain solve fails below this step size (throwing time-step crash
     error).

   INCLUDES:
   - ``[mpc-spec]`` *Is an* MPC_.

*/

#ifndef PKS_MPC_PERMAFROST_SPLIT_FLUX_HH_
#define PKS_MPC_PERMAFROST_SPLIT_FLUX_HH_

#include "PK.hh"
#include "mpc.hh"
#include "primary_variable_field_evaluator.hh"

namespace Amanzi {

class MPCPermafrostSplitFlux : public MPC<PK> {

 public:
  MPCPermafrostSplitFlux(Teuchos::ParameterList& FElist,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~MPCPermafrostSplitFlux() = default;

  // PK methods
  // -- initialize in reverse order
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;
  virtual void Setup(const Teuchos::Ptr<State>& S) override;

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;

  virtual double get_dt() override;
  virtual void set_dt(double dt) override;

  virtual void CommitStep(double t_old, double t_new,
                          const Teuchos::RCP<State>& S) override;
  virtual bool ValidStep() override;

 protected:
  bool AdvanceStep_Standard_(double t_old, double t_new, bool reinit);
  bool AdvanceStep_Subcycled_(double t_old, double t_new, bool reinit);

  void CopyPrimaryToStar_(const Teuchos::Ptr<const State>& S,
          const Teuchos::Ptr<State>& S_star);
  void CopyStarToPrimary_(double dt);

  void CopyPrimaryToStar_DomainSet_(const Teuchos::Ptr<const State>& S,
          const Teuchos::Ptr<State>& S_star);
  void CopyStarToPrimary_DomainSet_Pressure_(double dt);
  void CopyStarToPrimary_DomainSet_Flux_(double dt);
  void CopyStarToPrimary_DomainSet_Hybrid_(double dt);

  void CopyPrimaryToStar_Standard_(const Teuchos::Ptr<const State>& S,
          const Teuchos::Ptr<State>& S_star);
  void CopyStarToPrimary_Standard_Pressure_(double dt);
  void CopyStarToPrimary_Standard_Flux_(double dt);
  void CopyStarToPrimary_Standard_Hybrid_(double dt);


 protected:
  Key p_primary_variable_;
  Key p_primary_variable_suffix_;
  Key p_sub_primary_variable_;
  Key p_sub_primary_variable_suffix_;
  Key p_primary_variable_star_;
  Key p_conserved_variable_star_;
  Key p_lateral_flow_source_;
  Key p_lateral_flow_source_suffix_;

  Key T_primary_variable_;
  Key T_primary_variable_suffix_;
  Key T_sub_primary_variable_;
  Key T_sub_primary_variable_suffix_;
  Key T_primary_variable_star_;
  Key T_conserved_variable_star_;
  Key T_lateral_flow_source_;
  Key T_lateral_flow_source_suffix_;

  Key cv_key_;

  Key domain_set_;
  Key domain_;
  Key domain_sub_;
  Key domain_star_;
  Key domain_snow_;

  std::string coupling_;

  bool is_domain_set_;
  bool subcycled_;
  double subcycled_target_dt_;
  double subcycled_min_dt_;
  double cycle_dt_;

  // note, only one of these set will be used, the pointers or the vectors
  Teuchos::RCP<PrimaryVariableFieldEvaluator> p_eval_pvfe_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> T_eval_pvfe_;
  std::vector<Teuchos::RCP<PrimaryVariableFieldEvaluator>> p_eval_pvfes_;
  std::vector<Teuchos::RCP<PrimaryVariableFieldEvaluator>> T_eval_pvfes_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrostSplitFlux> reg_;
};

} // close namespace Amanzi

#endif
