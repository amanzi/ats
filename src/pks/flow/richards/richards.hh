/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! Two-phase, variable density Richards equation.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


/*!

Solves Richards equation:

.. math::
  \frac{\partial \Theta}{\partial t} - \nabla \frac{k_r n_l}{\mu} K ( \nabla p + \rho g \cdot \hat{z} ) = Q_w


Options:

Variable naming:

* `"domain`" ``[string]`` **""**  Defaults to the base subsurface domain.

* `"primary variable`" ``[string]`` The primary variable associated with this PK, typically `"pressure`"


Other variable names, typically not set as the default is basically always good:

* `"conserved quantity suffix`" ``[string]`` **"water_content"**  If set, changes the conserved quantity key.

* `"conserved quantity key`" ``[string]`` **"DOMAIN-CONSERVED_QUANTITY_SUFFIX"** Typically not set, default is good. ``[mol]``

* `"mass density key`" ``[string]`` **"DOMAIN-mass_density_liquid"** liquid water density ``[kg m^-3]``

* `"molar density key`" ``[string]`` **"DOMAIN-molar_density_liquid"** liquid water density ``[mol m^-3]``

* `"permeability key`" ``[string]`` **"DOMAIN-permeability"** permeability of the soil medium ``[m^2]``

* `"conductivity key`" ``[string]`` **"DOMAIN-relative_permeability"** scalar coefficient of the permeability ``[-]``

* `"upwind conductivity key`" ``[string]`` **"DOMAIN-upwind_relative_permeability"** upwinded (face-based) scalar coefficient of the permeability.  Note the units of this are strange, but this represents :math:`\frac{n_l k_r}{\mu}`  ``[mol kg^-1 s^1 m^-2]``

* `"darcy flux key`" ``[string]`` **"DOMAIN-mass_flux"** mass flux across a face ``[mol s^-1]``

* `"darcy flux direction key`" ``[string]`` **"DOMAIN-mass_flux_direction"** direction of the darcy flux (used in upwinding :math:`k_r`) ``[??]``

* `"darcy velocity key`" ``[string]`` **"DOMAIN-darcy_velocity"** darcy velocity vector, interpolated from faces to cells ``[m s^-1]``

* `"darcy flux key`" ``[string]`` **"DOMAIN-mass_flux"** mass flux across a face ``[mol s^-1]``

* `"saturation key`" ``[string]`` **"DOMAIN-saturation_liquid"** volume fraction of the liquid phase ``[-]``

Discretization control:

* `"diffusion`" ``[list]`` An OperatorDiffusion_ spec describing the (forward) diffusion operator

* `"diffusion preconditioner`" ``[list]`` An OperatorDiffusion_ spec describing the diffusive parts of the preconditioner.


Time integration and timestep control:

* `"initial time step`" ``[double]`` **1.** Max initial time step size ``[s]``

* `"time integrator`" ``[time-integrator-spec]`` is a TimeIntegrator_.  Note
  that this is only provided if this Richards PK is not strongly coupled to
  other PKs.

* `"linear solver`" ``[linear-solver-spec]`` is a LinearSolver_ spec.  Note
  that this is only used if this PK is not strongly coupled to other PKs.

* `"preconditioner`" ``[preconditioner-spec]`` is a Preconditioner_ spec.
  Note that this is only used if this PK is not strongly coupled to other PKs.

* `"initial condition`" ``[initial-condition-spec]`` See InitialConditions_.
  Additionally, the following parameter is supported:

  - `"initialize faces from cell`" ``[bool]`` **false** Indicates that the
    primary variable field has both CELL and FACE objects, and the FACE values
    are calculated as the average of the neighboring cells.

Error control:

* `"absolute error tolerance`" [double] **DERIVED** Defaults to a porosity of 0.5 * a saturation of 0.1 * n_l.  A small, but significant, amount of water.

* `"relative error tolerance`" [double] **1** Take the error relative to the amount of water present in that cell.

* `"flux tolerance`" [double] **1**

  Multiplies the error in flux (on a face) relative to the min of water in the
  neighboring cells.  Typically only changed if infiltration is very small and
  the boundary condition is not converging, at which point it can be decreased
  by an order of magnitude at a time until the boundary condition is
  satisfied.

Boundary conditions:

* `"boundary conditions`" ``[subsurface-flow-bc-spec]`` **defaults to Neuman, 0 normal flux**  See `Flow-specific Boundary Conditions`_

Physics control:

* `"permeability rescaling`" ``[double]`` **1** Typically 1e7 or order :math:`sqrt(K)` is about right.  This rescales things to stop from multiplying by small numbers (permeability) and then by large number (:math:`\rho / \mu`).

May inherit options from PKPhysicalBDFBase_

*/

#ifndef PK_FLOW_RICHARDS_HH_
#define PK_FLOW_RICHARDS_HH_

#include "wrm_partition.hh"
#include "BoundaryFunction.hh"
#include "upwinding.hh"

#include "LinearOperatorFactory.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorAccumulation.hh"

#include "PK_Factory.hh"
//#include "PK_PhysicalBDF_ATS.hh"
// #include "pk_factory_ats.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {

// forward declarations
class MPCSubsurface;
class PredictorDelegateBCFlux;
namespace WhetStone { class Tensor; }

namespace Flow {

class Richards : public PK_PhysicalBDF_Default {

public:

  Richards(Teuchos::ParameterList& FElist,
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~Richards() {}

  // virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {CalculateDiagnostics(S);};

  // -- Setup data.
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  // -- limit changes in a valid time step
  virtual bool ValidStep();

  // -- Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // problems with pressures -- setting a range of admissible pressures
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up);

  // evaluating consistent faces for given BCs and cell values
  virtual void CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u);

protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);
  virtual void SetupRichardsFlow_(const Teuchos::Ptr<State>& S);

  // boundary condition members
  virtual void UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S, bool kr=true);

  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual void SetAbsolutePermeabilityTensor_(const Teuchos::Ptr<State>& S);
  virtual bool UpdatePermeabilityData_(const Teuchos::Ptr<State>& S);
  virtual bool UpdatePermeabilityDerivativeData_(const Teuchos::Ptr<State>& S);

  virtual void UpdateVelocity_(const Teuchos::Ptr<State>& S);

  // physical methods
  // -- diffusion term
  virtual void ApplyDiffusion_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& g);

  // virtual void AddVaporDiffusionResidual_(const Teuchos::Ptr<State>& S,
  //         const Teuchos::Ptr<CompositeVector>& g);
  // virtual void ComputeVaporDiffusionCoef(const Teuchos::Ptr<State>& S, 
  //                                        Teuchos::RCP<CompositeVector>& vapor_diff, 
  //                                        std::string var_name);
 


  // -- accumulation term
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);

  // -- Add any source terms into the residual.
  virtual void AddSources_(const Teuchos::Ptr<State>& S,
                           const Teuchos::Ptr<CompositeVector>& f);
  virtual void AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h);
  
  // -- gravity contributions to matrix or vector
  // virtual void AddGravityFluxes_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
  //         const Teuchos::Ptr<const CompositeVector>& rel_perm,
  //         const Teuchos::Ptr<const CompositeVector>& rho,
  //         const Teuchos::Ptr<Operators::MatrixMFD>& matrix);

  // virtual void AddGravityFluxes_FV_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
  //         const Teuchos::Ptr<const CompositeVector>& rel_perm,
  //         const Teuchos::Ptr<const CompositeVector>& rho,
  //         const Teuchos::Ptr<Operators::Matrix_TPFA>& matrix);

  // virtual void AddGravityFluxesToVector_(const Teuchos::Ptr<const Epetra_Vector>& g_vec,
  //         const Teuchos::Ptr<const CompositeVector>& rel_perm,
  //         const Teuchos::Ptr<const CompositeVector>& rho,
  //         const Teuchos::Ptr<CompositeVector>& mass_flux);

  // Nonlinear version of CalculateConsistentFaces()
  // virtual void CalculateConsistentFacesForInfiltration_(
  //     const Teuchos::Ptr<CompositeVector>& u);
  virtual bool ModifyPredictorConsistentFaces_(double h, Teuchos::RCP<TreeVector> u);
  virtual bool ModifyPredictorWC_(double h, Teuchos::RCP<TreeVector> u);
  virtual bool ModifyPredictorFluxBCs_(double h, Teuchos::RCP<TreeVector> u);

  // virtual void PreconWC_(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Possibly modify the correction before it is applied
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> du);




  


  
protected:
  // control switches
  Operators::UpwindMethod Krel_method_;
  bool infiltrate_only_if_unfrozen_;
  bool modify_predictor_with_consistent_faces_;
  bool modify_predictor_wc_;
  bool symmetric_;
  bool precon_wc_;
  bool is_source_term_;
  bool source_term_is_differentiable_;
  bool explicit_source_;
  bool precon_used_;
  bool clobber_surf_kr_;
  bool clobber_boundary_flux_dir_;
  
  // coupling terms
  bool coupled_to_surface_via_head_; // surface-subsurface Dirichlet coupler
  bool coupled_to_surface_via_flux_; // surface-subsurface Neumann coupler

  // -- water coupler coupling parameters
  double surface_head_cutoff_;
  double surface_head_cutoff_alpha_;
  double surface_head_eps_;

  // permeability
  Teuchos::RCP<std::vector<WhetStone::Tensor> > K_;  // absolute permeability
  Teuchos::RCP<Operators::Upwinding> upwinding_;
  Teuchos::RCP<Operators::Upwinding> upwinding_deriv_;
  Teuchos::RCP<Flow::WRMPartition> wrms_;
  bool upwind_from_prev_flux_;

  // mathematical operators
  Teuchos::RCP<Operators::Operator> matrix_; // pc in PKPhysicalBDFBase
  Teuchos::RCP<Operators::OperatorDiffusionWithGravity> matrix_diff_;
  Teuchos::RCP<Operators::OperatorDiffusionWithGravity> preconditioner_diff_;
  Teuchos::RCP<Operators::OperatorDiffusionWithGravity> face_matrix_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> preconditioner_acc_;
  Teuchos::RCP<Operators::Operator> lin_solver_;

  // flag to do jacobian and therefore coef derivs
  bool jacobian_;
  int iter_;
  double iter_counter_time_;
  int jacobian_lag_;
  

  // residual vector for vapor diffusion
  Teuchos::RCP<CompositeVector> res_vapor;
  // note PC is in PKPhysicalBDFBase

  // custom enorm tolerances
  double flux_tol_;

  // boundary condition data
  Teuchos::RCP<Functions::BoundaryFunction> bc_pressure_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_head_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_seepage_infilt_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_infiltration_;

  // delegates
  bool modify_predictor_bc_flux_;
  bool modify_predictor_first_bc_flux_;
  Teuchos::RCP<PredictorDelegateBCFlux> flux_predictor_;

  // is this a dynamic mesh problem
  bool dynamic_mesh_;

  // is vapor turned on
  bool vapor_diffusion_;

  // scale for perm
  double perm_scale_;

  // limiters
  double p_limit_;
  double patm_limit_;

  // valid step controls
  double sat_change_limit_;
  double sat_ice_change_limit_;

  // keys
  Key mass_dens_key_;
  Key molar_dens_key_;
  Key perm_key_;
  Key coef_key_, dcoef_key_;
  Key uw_coef_key_, duw_coef_key_;
  Key flux_key_;
  Key flux_dir_key_;
  Key velocity_key_;
  Key source_key_;
  Key ss_flux_key_;
  Key sat_key_;
  Key sat_ice_key_;
  Key sat_gas_key_;

 private:
  // factory registration
  static RegisteredPKFactory<Richards> reg_;

  // Richards has a friend in couplers...
  friend class Amanzi::MPCSubsurface;

};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
