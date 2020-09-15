/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A coupler which solves flow and energy both surface and subsurface.


/*!

This MPC handles the coupling of surface energy and flow to subsurface energy
and flow for integrated hydrology with freeze/thaw processes.

.. _mpc-permafrost-spec:
.. admonition:: mpc-permafrost-spec

   * `"PKs order`" ``[Array(string)]`` The user supplies the names of the
     coupled PKs.  The order must be {subsurface_flow_pk, subsurface_energy_pk,
     surface_flow_pk, surface_energy_pk}.

   * `"subsurface domain name`" ``[string]`` **domain** 

   * `"surface domain name`" ``[string]`` **surface** 

   * `"mass exchange flux key`" ``[string]`` **SURFACE_DOMAIN-surface_subsurface_flux**

   * `"energy exchange flux key`" ``[string]`` **SURFACE_DOMAIN-surface_subsurface_energy_flux**

   * `"water delegate`" ``[mpc-delegate-water-spec]`` A `Coupled Water
     Globalization Delegate`_ spec.

   INCLUDES:

   - ``[mpc-subsurface-spec]`` *Is a* `Subsurface MPC`_
    
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

  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk FunctionalResidual().
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // -- Apply preconditioner to r and returns the result in Pr.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> r, Teuchos::RCP<TreeVector> Pr);

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> r,
                       Teuchos::RCP<const TreeVector> u, 
                       Teuchos::RCP<TreeVector> du);

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
  Teuchos::RCP<PrimaryVariableFieldEvaluator> mass_exchange_pvfe_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> energy_exchange_pvfe_;
  
  // off-diagonal terms
  Teuchos::RCP<Operators::PDE_Accumulation> dE_dp_surf_;

  
  // EWC delegate for the surface
  //  Teuchos::RCP<MPCDelegateEWC> surf_ewc_;

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

} // namespace


#endif
