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


.. todo:
    Nearly all variable name roots are hard-coded here, this should get updated.

*/


#ifndef PK_FLOW_OVERLAND_HEAD_HH_
#define PK_FLOW_OVERLAND_HEAD_HH_

#include "BoundaryFunction.hh"
#include "DynamicBoundaryFunction.hh"
#include "upwinding.hh"

#include "Operator.hh"
#include "PDE_Diffusion.hh"
#include "PDE_Accumulation.hh"

//#include "pk_factory_ats.hh"
//#include "pk_physical_bdf_base.hh"
#include "PK_Factory.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {

class MPCSurfaceSubsurfaceDirichletCoupler;

namespace Flow {

class OverlandConductivityModel;
class HeightModel;

//class OverlandPressureFlow : public PKPhysicalBDFBase {
class OverlandPressureFlow : public PK_PhysicalBDF_Default {
 public:
  OverlandPressureFlow(Teuchos::ParameterList& pk_tree,
                       const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~OverlandPressureFlow() {}

  // main methods
  // -- Initialize owned (dependent) variables.
  virtual void Setup() override;

  // -- Initialize owned (dependent) variables.
  virtual void Initialize() override;

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // -- Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Tag& tag) override;

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  void FunctionalResidual(double t_old,
                          double t_new,
                          Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> g) override;

  // applies preconditioner to u and returns the result in Pu
  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  virtual bool
  ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override;

  // evaluating consistent faces for given BCs and cell values
  virtual void CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u);

  // -- Possibly modify the correction before it is applied
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h,
                   Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override;

 protected:
  // setup methods
  virtual void SetupOverlandFlow_();
  virtual void SetupPhysicalEvaluators_();

  // boundary condition members
  void ComputeBoundaryConditions_(const Tag& tag);
  virtual void UpdateBoundaryConditions_(const Tag& tag);
  virtual void ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u,
                                        const Teuchos::Ptr<const CompositeVector>& elev);

  virtual void
  FixBCsForOperator_(const Tag& tag, const Teuchos::Ptr<Operators::PDE_Diffusion>& diff_op);
  virtual void FixBCsForPrecon_(const Tag& tag);

  // computational concerns in managing abs, rel perm
  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual bool UpdatePermeabilityDerivativeData_(const Tag& tag);
  virtual bool UpdatePermeabilityData_(const Tag& tag);

  // physical methods
  // -- diffusion term
  void ApplyDiffusion_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g);
  // -- accumulation term
  void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);
  // -- source terms
  void AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g);
  void AddSourcesToPrecon_(double h);

  void test_ApplyPreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  friend class Amanzi::MPCSurfaceSubsurfaceDirichletCoupler;

  enum FluxUpdateMode {
    UPDATE_FLUX_ITERATION = 0,
    UPDATE_FLUX_TIMESTEP = 1,
    UPDATE_FLUX_VIS = 2,
    UPDATE_FLUX_NEVER = 3
  };

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
  Key dcond_key_;
  Key duw_cond_key_;
  Key mass_dens_key_;
  Key molar_dens_key_;
  Key cv_key_;
  Key source_key_;
  Key source_molar_dens_key_;
  Key ss_flux_key_;

  // control switches
  bool standalone_mode_; // domain mesh == surface mesh
  FluxUpdateMode update_flux_;
  Operators::UpwindMethod upwind_method_;

  bool is_source_term_;
  bool source_term_is_differentiable_;
  bool source_only_if_unfrozen_;

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

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_zero_gradient_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_critical_depth_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_level_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_tidal_;
  Teuchos::RCP<Functions::DynamicBoundaryFunction> bc_dynamic_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_level_flux_lvl_, bc_level_flux_vel_;

  // needed physical models
  Teuchos::RCP<Flow::OverlandConductivityModel> cond_model_;

  // factory registration
  static RegisteredPKFactory<OverlandPressureFlow> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
