/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Overland flow using the diffusion wave equation.
/*!

Solves the diffusion wave equation for overland flow with pressure as a primary variable:

.. math::
  \frac{\partial \Theta}{\partial t} - \nabla n_l k \nabla h(p) = Q_w

Uses the PK type:
`"overland flow, pressure basis`"
  

.. _overland-pressure-spec:
.. admonition:: overland-pressure-spec

    Keys name variables:

    * `"domain`" ``[string]`` **"surface"**  Defaults to the extracted surface mesh.

    * `"primary variable`" ``[string]`` The primary variable associated with
      this PK, typically `"DOMAIN-pressure`" Note there is no default -- this
      must be provided by the user.

    * `"boundary conditions`" ``[list]`` Defaults to Neuman, 0 normal flux.

    * `"overland conductivity evaluator`" ``[list]``
      See `Overland Conductivity Evaluator`_.

    IF

    * `"source term`" ``[bool]`` **false** Is there a source term?

    THEN

    * `"source key`" ``[string]`` **DOMAIN-water_source** Typically
      not set, as the default is good. ``[m s^-1]`` or ``[mol s^-1]``
    * `"water source in meters`" ``[bool]`` **true** Is the source term in ``[m s^-1]``?
    * `"source term is differentiable`" ``[bool]`` **true** Can the source term
      be differentiated with respect to the primary variable?
    * `"explicit source term`" ``[bool]`` **false** Apply the source term from
      the previous time step.

    END

    Math and solver algorithm options:

    * `"diffusion`" ``[pde-diffusion-spec]`` The (forward) diffusion operator,
      see PDE_Diffusion_.

    * `"diffusion preconditioner`" ``[pde-diffusion-spec]`` **optional** The
      inverse of the diffusion operator.  See PDE_Diffusion_.  Typically this
      is only needed to set Jacobian options, as all others probably should
      match those in `"diffusion`", and default to those values.

    * `"absolute error tolerance`" ``[double]`` **550.** Defaults to 1 cm of
      water.  A small, but significant, amount of water.

    * `"limit correction to pressure change [Pa]`" ``[double]`` **-1** If > 0,
      this limits an iterate's max pressure change to this value.  Not usually
      helpful.

    * `"limit correction to pressure change when crossing atmospheric [Pa]`" ``[double]`` **-1**
      If > 0, this limits an iterate's max pressure change
      to this value when they cross atmospheric pressure.  Not usually helpful.

    * `"allow no negative ponded depths`" ``[bool]`` **false** Modifies all
      correction updates to ensure only positive ponded depth is allowed.

    * `"min ponded depth for velocity calculation`" ``[double]`` **1.e-2** For
      ponded depth below this height, declare the velocity 0.

    * `"min ponded depth for tidal bc`" ``[double]`` **0.02** Control on the
      tidal boundary condition.  TODO: This should live in the BC spec?

    INCLUDES:

    - ``[pk-physical-bdf-default-spec]`` A `PK: Physical and BDF`_ spec.

    Everything below this point is usually not provided by the user, but are
    documented here for completeness.

    Keys name variables:

    * `"conserved quantity key`" ``[string]`` **DOMAIN-water_content** Typically
      not set, as the default is good. ``[mol]``
    * `"elevation key`" ``[string]`` **DOMAIN-elevation** Typically
      not set, as the default is good. ``[mol]``
    * `"slope magnitude key`" ``[string]`` **DOMAIN-slope_magnitude** Typically
      not set, as the default is good. ``[mol]``

    Algorithmic parameters:

    * `"coupled to subsurface via flux`" ``[bool]`` **false** Set by MPC.
    * `"coupled to subsurface via head`" ``[bool]`` **false** Set by MPC.

    * `"accumulation preconditioner`" ``[pde-accumulation-spec]`` **optional**
      The inverse of the accumulation operator.  See PDE_Accumulation_.
      Typically not provided by users, as defaults are correct.

    EVALUATORS:

    - `"conserved quantity`"
    - `"water content`"
    - `"cell volume`"
    - `"surface_subsurface_flux`"
    - `"elevation`"
    - `"slope magnitude`"
    - `"overland_conductivity`"
    - `"ponded_depth`"
    - `"pres_elev`"
    - `"source`"

*/


#ifndef PK_FLOW_OVERLAND_HEAD_HH_
#define PK_FLOW_OVERLAND_HEAD_HH_

#include "upwinding.hh"

#include "PKFactory.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {

class MPCSurfaceSubsurfaceDirichletCoupler;

namespace Operators {
class PDE_Diffusion;
class PDE_Accumulation;
} // namespace Operators

namespace Flow {

class OverlandPressureFlow : public PK_PhysicalBDF_Default {
 public:
  OverlandPressureFlow(const Comm_ptr_type& comm,
                       Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& plist,
                       const Teuchos::RCP<State>& S);

  //
  // PK methods
  // ------------------------------------------------------------------
  // Parse the local parameter list and add entries to the global list
  virtual void ParseParameterList_() override;

  // Set requirements of data and evaluators
  virtual void Setup() override;

  // Initialize owned (dependent) variables.
  virtual void Initialize() override;

  // Finalize a step as successful at the given tag.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Tag& tag) override;

  // type info used in PKFactory
  static const std::string type;
  virtual const std::string& getType() const override { return type; }

  //
  // BDF1_TI methods
  // ------------------------------------------------------------------
  // Compute the non-linear functional g = g(t, u, du/dt)
  void FunctionalResidual(double t_old,
                          double t_new,
                          Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> g) override;

  // Apply the preconditioner to u and return the result in Pu
  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  // Updates the preconditioner by linearizing at up
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // Globalization -- modify the predictor u to make a better guess for the new time
  virtual bool
  ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override;

  // Possibly modify the correction du before it is applied
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h,
                   Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override;

 protected:

  //
  // Protected, internal methods for better granularity of design
  // ------------------------------------------------------------------
  virtual void SetupOverlandFlow_();
  virtual void SetupPhysicalEvaluators_();

  // boundary condition members
  virtual void UpdateBoundaryConditions_(const Tag& tag);

  virtual void
  FixBCsForOperator_(const Tag& tag, const Teuchos::Ptr<Operators::PDE_Diffusion>& diff_op);
  virtual void FixBCsForPrecon_(const Tag& tag);

  // -- builds the upwinded conductivity
  virtual bool UpdatePermeabilityData_(const Tag& tag);
  virtual bool UpdatePermeabilityDerivativeData_(const Tag& tag);
  void ApplyDirichletBCs_(const Operators::BCs& bcs, CompositeVector& u, const CompositeVector& elev);

  // physical methods
  // -- diffusion term
  void ApplyDiffusion_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g);
  // -- accumulation term
  void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);
  // -- source terms
  void AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g);

  void test_ApplyPreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  friend class Amanzi::MPCSurfaceSubsurfaceDirichletCoupler;

  // keys
  Key potential_key_;
  Key flux_key_;
  Key flux_dir_key_;
  Key velocity_key_;
  Key elev_key_;
  Key pd_key_;
  Key pd_bar_key_;
  Key wc_bar_key_;
  Key cond_key_;
  Key uw_cond_key_;
  Key duw_cond_key_;
  Key mass_dens_key_;
  Key molar_dens_key_;
  Key source_key_;
  Key source_molar_dens_key_;
  Key ss_flux_key_;

  // control switches
  Operators::UpwindMethod upwind_method_;

  bool is_source_term_;
  bool source_in_meters_;
  // bool source_only_if_unfrozen_;

  bool modify_predictor_with_consistent_faces_;
  bool symmetric_;
  bool perm_update_required_;

  double p_limit_;
  double patm_limit_;
  bool patm_hard_limit_;
  double min_vel_ponded_depth_, min_tidal_bc_ponded_depth_;

  // coupling term
  bool coupled_to_subsurface_via_head_;
  bool coupled_to_subsurface_via_flux_;

  // newton correction
  bool jacobian_;
  int iter_;
  double iter_counter_time_;
  int jacobian_lag_;

  // work data space
  Teuchos::RCP<Operators::Upwinding> upwinding_;
  Teuchos::RCP<Operators::Upwinding> upwinding_dkdp_;

  // mathematical operators
  Teuchos::RCP<Operators::Operator> matrix_; // pc in PKPhysicalBDFBase
  Teuchos::RCP<Operators::PDE_Diffusion> matrix_diff_;
  Teuchos::RCP<Operators::PDE_Diffusion> face_matrix_diff_;
  Teuchos::RCP<Operators::PDE_Diffusion> preconditioner_diff_;
  Teuchos::RCP<Operators::PDE_Accumulation> preconditioner_acc_;

  bool precon_used_;
  bool precon_scaled_;

  // factory registration
  static RegisteredPKFactory<OverlandPressureFlow> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
