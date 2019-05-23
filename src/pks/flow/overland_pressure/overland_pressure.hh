/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! Overland flow using the diffusion wave equation.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

Solves the diffusion wave equation for overland flow with pressure as a primary variable:

.. math::
  \frac{\partial \Theta}{\partial t} - \nabla n_l k \nabla h(p) = Q_w


Options:

Variable naming:

* `"domain`" ``[string]`` **"surface"**  Defaults to the extracted surface mesh.

* `"primary variable`" ``[string]`` The primary variable associated with this PK, typically `"DOMAIN-pressure`"


Other variable names, typically not set as the default is basically always good:

* `"conserved quantity key`" ``[string]`` **"DOMAIN-water_content"** The conserved quantity name.

Discretization control:

* `"diffusion`" ``[list]`` An PDE_Diffusion_ spec describing the (forward) diffusion operator

* `"diffusion preconditioner`" ``[list]`` An PDE_Diffusion_ spec describing the diffusive parts of the preconditioner.

Time integration and timestep control:

* `"initial time step`" ``[double]`` **1.** Max initial time step size ``[s]``.

* `"time integrator`" ``[implicit-time-integrator-typed-spec]`` **optional** is a TimeIntegrator_ spec.
  Note that this is only used if this PK is not strongly coupled to other PKs.

* `"linear solver`" ``[linear-solver-typed-spec]`` **optional** is a LinearSolver_ spec.  Note
  that this is only used if this PK is not strongly coupled to other PKs.

* `"preconditioner`" ``[preconditioner-typed-spec]`` **optional** is a Preconditioner_ spec.
  Note that this is only used if this PK is not strongly coupled to other PKs.

* `"initial condition`" ``[initial-conditions-spec]`` See InitialConditions_.

Error control:

* `"absolute error tolerance`" ``[double]`` **550.** Defaults to 1 cm of water.  A small, but significant, amount of water.

* `"relative error tolerance`" ``[double]`` **1** Take the error relative to the amount of water present in that cell.

* `"flux tolerance`" ``[double]`` **1** Multiplies the error in flux (on a face)
  relative to the min of water in the neighboring cells.  Typically only
  changed if infiltration is very small and the boundary condition is not
  converging, at which point it can be decreased by an order of magnitude at a
  time until the boundary condition is satisfied.

Boundary conditions:

xx* `"boundary conditions`" ``[surface-flow-bc-spec]`` Defaults to Neuman, 0 normal flux.


May inherit options from PKPhysicalBDFBase_.

*/


#ifndef PK_FLOW_OVERLAND_HEAD_HH_
#define PK_FLOW_OVERLAND_HEAD_HH_

#include "BoundaryFunction.hh"
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
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  // -- Update diagnostics for vis.
  //virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);


  // evaluating consistent faces for given BCs and cell values
  virtual void CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u);
  
  // -- Compute a norm on u-du and return the result.
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du);

  // -- Possibly modify the correction before it is applied
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> du);

protected:
  // setup methods
  virtual void SetupOverlandFlow_(const Teuchos::Ptr<State>& S);
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

  // boundary condition members
  virtual void UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S);
  virtual void ApplyBoundaryConditions_(const Teuchos::Ptr<CompositeVector>& u,
          const Teuchos::Ptr<const CompositeVector>& elev);

  virtual void FixBCsForOperator_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<Operators::PDE_Diffusion>& diff_op);
  //virtual void FixBCsForFluxDirection_(const Teuchos::Ptr<State>& S);
  virtual void FixBCsForPrecon_(const Teuchos::Ptr<State>& S);
  // virtual void FixBCsForConsistentFaces_(const Teuchos::Ptr<State>& S);

  // computational concerns in managing abs, rel perm
  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual bool UpdatePermeabilityDerivativeData_(const Teuchos::Ptr<State>& S);
  virtual bool UpdatePermeabilityData_(const Teuchos::Ptr<State>& S);

  // physical methods
  // -- diffusion term
  void ApplyDiffusion_(const Teuchos::Ptr<State>& S,const Teuchos::Ptr<CompositeVector>& g);
  // -- accumulation term
  void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);
  // -- source terms
  void AddSourceTerms_(const Teuchos::Ptr<CompositeVector>& g);

  void test_ApplyPreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  friend class Amanzi::MPCSurfaceSubsurfaceDirichletCoupler;

  enum FluxUpdateMode {
    UPDATE_FLUX_ITERATION = 0,
    UPDATE_FLUX_TIMESTEP = 1,
    UPDATE_FLUX_VIS = 2,
    UPDATE_FLUX_NEVER = 3
  };

  // control switches
  bool standalone_mode_; // domain mesh == surface mesh
  FluxUpdateMode update_flux_;
  Operators::UpwindMethod upwind_method_;

  bool is_source_term_;
  Key source_key_;
  bool source_in_meters_;
  bool source_only_if_unfrozen_;
  
  bool modify_predictor_with_consistent_faces_;
  bool symmetric_;
  bool perm_update_required_;
  
  double p_limit_;
  double patm_limit_;
  bool patm_hard_limit_;

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
  Teuchos::RCP<Operators::Operator> lin_solver_;

  bool precon_used_;

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_zero_gradient_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_critical_depth_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_level_;
  
  // needed physical models
  Teuchos::RCP<Flow::OverlandConductivityModel> cond_model_;

  int niter_;

  // factory registration
  static RegisteredPKFactory<OverlandPressureFlow> reg_;
};

}  // namespace Flow
}  // namespace Amanzi

#endif
