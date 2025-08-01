/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

This MPC solves the thermal equivalent of :ref:`Operator Split Integrated
Hydrology`.  This is the permafrost analog, so deals with energy as well in a
similar strategy.  In this case advection and diffusion of energy are handled
in the "star" system, then the "primary" system solves the :ref:`Integrated
Thermal Hydrology` system without lateral fluxes and with extra coupling terms
equivalent to those for flow from :ref:`Operaor Split Integrated Hydrology`.

Like the flow-only case, this can be used with either a 3D subsurface solve, by
setting the 2nd sub-PK to be a 3D permafrost MPC, or a bunch of columns, but
setting the 2nd sub-PK to be a DomainSetMPC.

This is the model from `Jan et al Comp. Geosci. 2018
<https://doi.org/10.1007/s10596-017-9679-3>`_

`"PK type`" = `"operator split permafrost`"

.. _pk-operator-split-permafrost:
.. admonition:: pk-operator-split-permafrost

   * `"PKs order`" ``[Array(string)]`` Order is {star_system, primary_system}.
     Note that the sub-PKs are likely a :ref:`Surface Flow & Energy` MPC for
     the star system and a :ref:`Integrated Thermal Hydrology` MPC for the
     primary system.  Because the latter does not include lateral fluxes in the
     surface, the sub-PKs of the primary system may be :ref:`General Surface
     Balance` equations rather than :ref:`Overland Flow PK` and :ref:`Overland
     Energy with Ice` (thereby removing the lateral flow and energy transport
     components).

   * `"domain name`" ``[string]`` The subsurface domain, e.g. "domain" (for a
     3D subsurface ) or "column:*" (for the intermediate scale model.

   * `"star domain name`" ``[string]`` The surface domain, typically
     `"surface_star`" by convention.

   * `"coupling type`" ``[string]`` **hybrid** One of: `"pressure`" (pass the
     pressure field when coupling flow in the operator splitting), `"flux`"
     (pass the divergence of fluxes as a source), or `"hybrid`" a mixture of
     the two that seems the most robust.

   INCLUDES:
   - ``[pk-mpc-spec]`` *Is an* :ref:`MPC`.
   - ``[pk-subcycling-mpc-spec]`` *Is a* :ref:`Subcycling MPC`

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
  // -- read input plist
  virtual void parseParameterList() override;

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

  Tag get_ds_tag_next_(const std::string& subdomain);
  Tag get_ds_tag_current_(const std::string& subdomain);

 protected:
  Key p_primary_variable_;
  Key p_primary_variable_suffix_;
  Key p_sub_primary_variable_;
  Key p_sub_primary_variable_suffix_;
  Key p_primary_variable_star_;
  Key p_conserved_variable_;
  Key p_conserved_variable_suffix_;
  Key p_conserved_variable_star_;
  Key p_lateral_flow_source_;
  Key p_lateral_flow_source_suffix_;

  Key T_primary_variable_;
  Key T_primary_variable_suffix_;
  Key T_sub_primary_variable_;
  Key T_sub_primary_variable_suffix_;
  Key T_primary_variable_star_;
  Key T_conserved_variable_;
  Key T_conserved_variable_suffix_;
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
  bool ds_is_subcycling_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrostSplitFlux> reg_;
};

} // namespace Amanzi
