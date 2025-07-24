/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*!

This MPC handles the coupling of surface energy and flow to subsurface energy
and flow for integrated thermal hydrology with freeze/thaw processes.

Coupling for flow is identical to that in :ref:`Integrated Hydrology`.
Coupling between surface flow and energy is identical to that of :ref:`Surface
Flow & Energy`, while coupling between subsurface flow and energy is identical
to that of :ref:`Subsurface Flow & Energy`.  Coupling of energy between the
surface and subsurface is similar to flow; both continuity of temperature and
the diffusive flux of energy are enforced, while advected fluxes of energy are
passed with mass (exfiltration) flux.

This is the implementation of the model described in `Painter et al WRR 2016
<https://doi.org/10.1002/2015WR018427>`_.

`"PK type`" = `"permafrost model`"

.. _pk-permafrost-model-spec:
.. admonition:: pk-permafrost-model-spec

   * `"PKs order`" ``[Array(string)]`` Order must be {subsurface_flow_pk,
     subsurface_energy_pk, surface_flow_pk, surface_energy_pk}.

   * `"subsurface domain name`" ``[string]`` **domain**

   * `"surface domain name`" ``[string]`` **surface**

   * `"mass exchange flux key`" ``[string]`` **SURFACE_DOMAIN-surface_subsurface_flux**

   * `"energy exchange flux key`" ``[string]`` **SURFACE_DOMAIN-surface_subsurface_energy_flux**

   * `"supress Jacobian terms: d div surface q / dT`" ``[bool]`` **false** Turn
     off the Jacobian term associated with the derivative of overland
     conductivity with respect to temperature.

   * `"water delegate`" ``[mpc-delegate-water-spec]`` A :ref:`Coupled Water
     Globalization Delegate` spec.

   * `"surface ewc delegate`" ``[mpc-delegate-ewc-spec]`` A :ref:`EWC
     Globalization Delegate` spec.

   INCLUDES:

   - ``[mpc-subsurface-spec]`` *Is a* :ref:`Subsurface Flow & Energy`

 */

#ifndef PKS_MPC_PERMAFROST_FOUR_HH_
#define PKS_MPC_PERMAFROST_FOUR_HH_

#include "mpc_delegate_ewc.hh"
#include "mpc_delegate_water.hh"
#include "mpc_subsurface.hh"

namespace Amanzi {

class MPCPermafrost : public MPCSubsurface {
 public:
  MPCPermafrost(Teuchos::ParameterList& FElist,
                const Teuchos::RCP<Teuchos::ParameterList>& plist,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& soln);

  // call to allow a PK to modify its own list or lists of its children.
  virtual void parseParameterList() override;

  virtual void set_tags(const Tag& tag_current, const Tag& tag_next) override;

  virtual void Setup() override;
  virtual void Initialize() override;
  //virtual void set_tags(const Tag& tag_current, const Tag& tag_next);

  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk FunctionalResidual().
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g) override;

  // -- Apply preconditioner to r and returns the result in Pr.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> r,
                                  Teuchos::RCP<TreeVector> Pr) override;

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h,
                               Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) override;

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult ModifyCorrection(
    double h,
    Teuchos::RCP<const TreeVector> r,
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> du) override;

 protected:
  // sub PKs
  Teuchos::RCP<PK_PhysicalBDF_Default> domain_flow_pk_;
  Teuchos::RCP<PK_PhysicalBDF_Default> domain_energy_pk_;
  Teuchos::RCP<PK_PhysicalBDF_Default> surf_flow_pk_;
  Teuchos::RCP<PK_PhysicalBDF_Default> surf_energy_pk_;

  // sub meshes
  Key domain_surf_;
  Key domain_subsurf_;
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;

  // Primary variable evaluators for exchange fluxes
  Key mass_exchange_key_;
  Key energy_exchange_key_;

  // off-diagonal terms
  // -- d ( dE/dt ) / dp terms
  Teuchos::RCP<Operators::PDE_Accumulation> dE_dp_surf_;
  // -- d ( div q ) / dT  terms
  Teuchos::RCP<Operators::PDE_Diffusion> ddivq_dT_;

  Key surf_temp_key_;
  Key surf_pres_key_;
  Key surf_e_key_;
  Key surf_wc_key_;
  Key surf_tc_key_;
  Key surf_kr_key_;
  Key surf_kr_uw_key_;
  Key surf_potential_key_;
  Key surf_pd_key_;
  Key surf_enth_key_;
  Key surf_water_flux_key_;

  // EWC delegate for the surface
  Teuchos::RCP<MPCDelegateEWC> surf_ewc_;

  // Water delegate
  Teuchos::RCP<MPCDelegateWater> water_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCPermafrost> reg_;

  Key domain_surf, domain_ss;
};

} // namespace Amanzi


#endif
