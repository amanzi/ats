/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A coupler which integrates surface and subsurface flow implicitly, without
//! relying on injecting the surface operator into the subsurface operator.
/*!

Couples Richards equation to surface water through continuity of both pressure
and fluxes, as does MPCCoupledWater, but WITHOUT relying on MFD's face-based
structure and the coincidence of surface cells with subsurface top boundary
faces. This works with any subsurface discretization, not just MFD's face+cell
discretization. The subsurface problem is solved with the surface pressure
imposed directly as a Dirichlet BC at the top face ("coupled to surface via
head" on the subsurface PK), the resulting flux becomes a source for the
surface equation ("coupled to subsurface via head" on the surface PK), and the
preconditioner is a genuine 2x2 TreeOperator with explicit off-diagonal
cross-mesh coupling blocks.

.. _mpc-coupled-water-noninjected-spec:
.. admonition:: mpc-coupled-water-noninjected-spec

   * `"PKs order`" ``[Array(string)]`` The user supplies the names of the
     coupled PKs.  The order must be {subsurface_flow_pk, surface_flow_pk}
     (subsurface first).

   * `"subsurface domain name`" ``[string]`` **domain**

   * `"surface domain name`" ``[string]`` **surface**

   * `"exfiltration key`" ``[string]`` **SURFACE_DOMAIN-surface_subsurface_flux**
     Key of a diagnostic field (always signed so that positive is
     exfiltration, subsurface to surface) storing the exchange flux computed
     each residual evaluation.

   * `"subsurface darcy flux key`" ``[string]`` **DOMAIN-water_flux**
   * `"subsurface upwind coefficient key`" ``[string]`` **DOMAIN-upwind_conductivity**

   * `"include Newton correction`" ``[bool]`` **false** Single switch for all
     dkr/dp Jacobian terms: Richards' own diagonal term (forwarded into its
     own "diffusion preconditioner" sublist) and this MPC's off-diagonal
     coupling-block cross terms. Always either both on or both off.

   INCLUDES:

   - ``[strong-mpc-spec]`` *Is a* `Strong MPC`_

*/

#ifndef PKS_MPC_COUPLED_WATER_NONINJECTED_HH_
#define PKS_MPC_COUPLED_WATER_NONINJECTED_HH_

#include "Operator.hh"
#include "Operator_Schema.hh"
#include "Op_MeshInjection.hh"
#include "Op_Cell_Cell.hh"
#include "TreeOperator.hh"
#include "PDE_DiffusionWithGravity.hh"

#include "PK_PhysicalBDF_Default.hh"
#include "MPCStrong.hh"
#include "mpc_delegate_water.hh"

namespace Amanzi {

class MPCCoupledWaterNonInjected : public MPCStrong<PK_PhysicalBDF_Default> {
 public:
  MPCCoupledWaterNonInjected(const Comm_ptr_type& comm,
                            Teuchos::ParameterList& pk_tree,
                            const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                            const Teuchos::RCP<State>& S);

  // Parse the local parameter list and add entries to the global list
  virtual void parseParameterList() override;

  virtual void setup() override;

  // Initialize: ensure continuity of ICs between subsurface and surface
  // pressure before the first residual/preconditioner evaluation (mirrors
  // MPCCoupledWater::initialize()).
  virtual void initialize() override;

  // type info used in PKFactory
  virtual const std::string& getType() const override { return pk_type_; }

  // -- computes the non-linear functional g = g(t,u,udot)
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g) override;

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  // -- Modify the predictor.
  virtual bool
  ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override;

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h,
                   Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override;

 protected:
  // Assembles the (0,1)/(1,0) off-diagonal coupling blocks of precon_. Also
  // folds the (0,0) surface-source Jacobian addition (d(exchange
  // flux)/d(p_surf)) into block (1,1)'s existing diagonal
  void UpdateCouplingBlocks_();

 protected:
  static const std::string pk_type_;

  std::string domain_surf_, domain_ss_;
  Key exfilt_key_;
  Key ss_uw_coef_key_;
  Key ss_flux_key_;
  Key ss_pressure_key_;

  // Single on/off switch for all Newton-correction dkr/dp terms: Richards'
  // own diagonal term (forwarded into its "diffusion preconditioner"
  // sublist as "include Newton correction") AND the coupling blocks'
  // cross terms computed in UpdateCouplingBlocks_() below. Set from this
  // MPC's own plist in parseParameterList(), so the two can never disagree.
  bool include_newton_correction_;

  // sub PKs
  Teuchos::RCP<PK_PhysicalBDF_Default> domain_flow_pk_;
  Teuchos::RCP<PK_PhysicalBDF_Default> surf_flow_pk_;

  // sub meshes
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;

  // coupled preconditioner: a 2x2 TreeOperator, with the two sub-PKs' own
  // preconditioners on the diagonal and rectangular, cross-mesh off-diagonal
  // coupling blocks
  Teuchos::RCP<Operators::TreeOperator> precon_;

  // off-diagonal coupling blocks: (0,1) = d(subsurface residual)/d(p_surf),
  // (1,0) = d(surface residual)/d(p_sub). Both share a single injection map
  // (surface-cell-local-index -> subsurface-top-cell GID), built once in
  // setup(); only their `diag` scalar coefficients are updated each
  // UpdateCouplingBlocks_() call.
  Teuchos::RCP<Operators::Operator_Schema> coupling_op_01_, coupling_op_10_;
  Teuchos::RCP<Operators::Op_MeshInjection> coupling_local_op_01_, coupling_local_op_10_;

  // (1,1) diagonal addition for overland's own source-Jacobian gap
  // (d(surface residual)/d(p_surf) via the exchange flux)
  Teuchos::RCP<Operators::Op_Cell_Cell> coupling_acc_surf_;

  // water globalization
  Teuchos::RCP<MPCDelegateWater> water_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;

  // Set true after precon_->SymbolicAssembleMatrix() has been called once
  // (lazily, on the first UpdatePreconditioner() call) -- see there.
  bool symbolic_assembled_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCCoupledWaterNonInjected> reg_;
};

} // namespace Amanzi


#endif
