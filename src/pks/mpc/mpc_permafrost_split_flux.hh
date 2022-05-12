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

   INCLUDES:
   - ``[mpc-spec]`` *Is an* MPC_.
   - ``[mpc-subcycled-spec]`` *Is a* MPCSubcycled_

*/

#pragma once

#include "PK.hh"
#include "mpc_subcycled.hh"

namespace Amanzi {

class MPCPermafrostSplitFlux : public MPCSubcycled {

 public:
  MPCPermafrostSplitFlux(Teuchos::ParameterList& FElist,
          const Teuchos::RCP<Teuchos::ParameterList>& plist,
          const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& solution);

  // PK methods
  // -- initialize in reverse order
  virtual void Initialize() override;
  virtual void Setup() override;

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

 protected:

  void CopyPrimaryToStar_(const Tag& prim, const Tag& star);
  void CopyStarToPrimary_(const Tag& star_current, const Tag& star_next,
          const Tag& prim_current, const Tag& prim_next);

  void CopyPrimaryToStar_DomainSet_(const Tag& prim, const Tag& star);
  void CopyStarToPrimary_DomainSet_Pressure_(const Tag& star_current, const Tag& star_next,
          const Tag& prim_current, const Tag& prim_next);
  void CopyStarToPrimary_DomainSet_Flux_(const Tag& star_current, const Tag& star_next,
          const Tag& prim_current, const Tag& prim_next);
  void CopyStarToPrimary_DomainSet_Hybrid_(const Tag& star_current, const Tag& star_next,
          const Tag& prim_current, const Tag& prim_next);

  void CopyPrimaryToStar_Standard_(const Tag& prim, const Tag& star);
  void CopyStarToPrimary_Standard_Pressure_(const Tag& star_current, const Tag& star_next,
          const Tag& prim_current, const Tag& prim_next);
  void CopyStarToPrimary_Standard_Flux_(const Tag& star_current, const Tag& star_next,
          const Tag& prim_current, const Tag& prim_next);
  void CopyStarToPrimary_Standard_Hybrid_(const Tag& star_current, const Tag& star_next,
          const Tag& prim_current, const Tag& prim_next);


 protected:
  Key p_primary_variable_;
  Key p_primary_variable_suffix_;
  Key p_sub_primary_variable_;
  Key p_sub_primary_variable_suffix_;
  Key p_primary_variable_star_;
  Key p_conserved_variable_;
  Key p_conserved_variable_star_;
  Key p_lateral_flow_source_;
  Key p_lateral_flow_source_suffix_;

  Key T_primary_variable_;
  Key T_primary_variable_suffix_;
  Key T_sub_primary_variable_;
  Key T_sub_primary_variable_suffix_;
  Key T_primary_variable_star_;
  Key T_conserved_variable_;
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

 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrostSplitFlux> reg_;
};

} // close namespace Amanzi


