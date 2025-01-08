/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! An operator-split coupled_water coupler, splitting overland flow from subsurface.
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

Note that this can be used with either a 3D subsurface solve, by setting the
2nd sub-PK to be a 3D permafrost MPC, or a bunch of columns, but setting the
2nd sub-PK to be a DomainSetMPC.


.. _mpc_permafrost_split_flux-spec
.. admonition:: mpc_permafrost_split_flux-spec

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
   - ``[mpc_subcycled-spec]`` *Is a* MPCSubcycled_

*/

#pragma once

#include "PK.hh"
#include "mpc_subcycled.hh"

namespace Amanzi {

class MPCCoupledWaterSplitFlux : public MPCSubcycled {
 public:
  MPCCoupledWaterSplitFlux(Teuchos::ParameterList& FElist,
                           const Teuchos::RCP<Teuchos::ParameterList>& plist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution);

  // PK methods
  void parseParameterList() override;

  // -- initialize in reverse order
  virtual void Initialize() override;
  virtual void Setup() override;

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

 protected:
  void CopyPrimaryToStar_();
  void CopyStarToPrimary_();

  void CopyPrimaryToStar_DomainSet_();
  void CopyStarToPrimary_DomainSet_Pressure_();
  void CopyStarToPrimary_DomainSet_Flux_();
  void CopyStarToPrimary_DomainSet_Hybrid_();

  void CopyPrimaryToStar_Standard_();
  void CopyStarToPrimary_Standard_Pressure_();
  void CopyStarToPrimary_Standard_Flux_();
  void CopyStarToPrimary_Standard_Hybrid_();

  Tag get_ds_tag_next_(const std::string& subdomain)
  {
    if (subcycling_[1])
      return Tag{ Keys::cleanName(tags_[1].second.get() + "_" +
                                  Keys::getDomainSetIndex(subdomain)) };
    else
      return Tag{ tags_[1].second };
  }
  Tag get_ds_tag_current_(const std::string& subdomain)
  {
    if (subcycling_[1])
      return Tag{ Keys::cleanName(tags_[1].first.get() + "_" +
                                  Keys::getDomainSetIndex(subdomain)) };
    else
      return Tag{ tags_[1].first };
  }


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

  Key cv_key_;

  Key domain_set_;
  Key domain_;
  Key domain_sub_;
  Key domain_star_;

  std::string coupling_;

  bool is_domain_set_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCCoupledWaterSplitFlux> reg_;
};

} // namespace Amanzi
