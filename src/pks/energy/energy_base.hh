/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! An advection-diffusion equation for energy.

/*!

Solves an advection-diffusion equation for energy:

.. math::
    \frac{\partial E}{\partial t} - \nabla \cdot \kappa \nabla T + \nabla \cdot \mathbf{q} e(T) = Q_w e(T) + Q_e

.. todo:: Document the energy error norm!
    
.. _energy-pk-spec:
.. admonition:: energy-pk-spec

    * `"domain`" ``[string]`` **"domain"**  Defaults to the subsurface mesh.

    * `"primary variable`" ``[string]`` The primary variable associated with
      this PK, typically `"DOMAIN-temperature`" Note there is no default -- this
      must be provided by the user.

    * `"boundary conditions`" ``[energy-bc-spec]`` Defaults to 0 diffusive flux
      boundary condition.  See `Energy-specific Boundary Conditions`_
      
    * `"thermal conductivity evaluator`"
      ``[thermal-conductivity-evaluator-spec]`` The thermal conductivity.  This
      needs to go away, and should get moved to State.

    * `"absolute error tolerance`" ``[double]`` **76.e-6** A small amount of
      energy, see error norm. `[MJ]`
      
    * `"upwind conductivity method`" ``[string]`` **arithmetic mean** Method of
      moving cell-based thermal conductivities onto faces.  One of:

      - `"arithmetic mean`" the default, average of neighboring cells
      - `"cell centered`" harmonic mean

    IF
      
    * `"explicit advection`" ``[bool]`` **false** Treat the advection term implicitly.
    
    ELSE

    * `"supress advective terms in preconditioner`" ``[bool]`` **false**
      Typically subsurface energy equations are strongly diffusion dominated,
      and the advective terms may add little.  With this flag on, we ignore
      theem in the preconditioner, making an easier linear solve and often not
      negatively impacting the nonlinear solve.

    * `"advection preconditioner`" ``[pde-advection-spec]`` **optional**
      Typically defaults are correct.

    END
      
    * `"diffusion`" ``[pde-diffusion-spec]`` See PDE_Diffusion_, the diffusion operator.

    * `"diffusion preconditioner`" ``[pde-diffusion-spec]`` See
      PDE_Diffusion_, the inverse operator.  Typically only adds Jacobian
      terms, as all the rest default to those values from `"diffusion`".

    IF
    
    * `"source term`" ``[bool]`` **false** Is there a source term?
    
    THEN
    
    * `"source key`" ``[string]`` **DOMAIN-total_energy_source** Typically
      not set, as the default is good. ``[MJ s^-1]``

    * `"source term is differentiable`" ``[bool]`` **true** Can the source term
      be differentiated with respect to the primary variable?

    * `"source term finite difference`" ``[bool]`` **false** If the source term
      is not diffferentiable, we can do a finite difference approximation of
      this derivative anyway.  This is useful for difficult-to-differentiate
      terms like a surface energy balance, which includes many terms.

    EVALUATORS:

    - `"source term`"
    
    END

    Globalization:
    
    * `"modify predictor with consistent faces`" ``[bool]`` **false** In a
      face+cell diffusion discretization, this modifies the predictor to make
      sure that faces, which are a DAE, are consistent with the predicted cells
      (i.e. face fluxes from each sides match).

    * `"modify predictor for freezing`" ``[bool]`` **false** A simple limiter
      that keeps temperature corrections from jumping over the phase change.

    * `"limit correction to temperature change [K]`" ``[double]`` **-1.0** If >
      0, stops nonlinear updates from being too big through clipping.

    The following are rarely set by the user, as the defaults are typically right.
    
    Variable names:

    * `"conserved quantity key`" ``[string]`` **DOMAIN-energy** The total energy :math:`E` `[MJ]`
    * `"energy key`" ``[string]`` **DOMAIN-energy** The total energy :math:`E`, also the conserved quantity. `[MJ]`
    * `"water content key`" ``[string]`` **DOMAIN-water_content** The total mass :math:`\Theta`, used in error norm `[mol]`
    * `"enthalpy key`" ``[string]`` **DOMAIN-enthalpy** The specific enthalpy :math`e` `[MJ mol^-1]`
    * `"flux key`" ``[string]`` **DOMAIN-mass_flux** The mass flux :math:`\mathbf{q}` used in advection. `[mol s^-1]`
    * `"diffusive energy flux`" ``[string]`` **DOMAIN-diffusive_energy_flux** :math:`\mathbf{q_e}` `[MJ s^-1]`
    * `"advected energy flux`" ``[string]`` **DOMAIN-advected_energy_flux** :math:`\mathbf{q_e^{adv}} = q e` `[MJ s^-1]`
    * `"thermal conductivity`" ``[string]`` **DOMAIN-thermal_conductivity** Thermal conductivity on cells `[W m^-1 K^-1]`
    * `"upwinded thermal conductivity`" ``[string]`` **DOMAIN-upwinded_thermal_conductivity** Thermal conductivity on faces `[W m^-1 K^-1]`
    
    * `"advection`" ``[pde-advection-spec]`` **optional** The PDE_Advection_ spec.  Only one current implementation, so defaults are typically fine.
    
    * `"accumulation preconditioner`" ``[pde-accumulation-spec]`` **optional**
      The inverse of the accumulation operator.  See PDE_Accumulation_.
      Typically not provided by users, as defaults are correct.

    IF
    
    * `"coupled to surface via flux`" ``[bool]`` **false** If true, apply
      surface boundary conditions from an exchange flux.  Note, if this is a
      coupled problem, it is probably set by the MPC.  No need for a user to
      set it.
      
    THEN

    * `"surface-subsurface energy flux key`" ``[string]`` **DOMAIN-surface_subsurface_energy_flux**
    
    END

    * `"coupled to surface via temperature`" ``[bool]`` **false** If true, apply
      surface boundary conditions from the surface temperature (Dirichlet).


    EVALUATORS:

    - `"enthalpy`"
    - `"cell volume`"
    - `"thermal conductivity`"
    - `"conserved quantity`"
    - `"energy`"
    
 */


#ifndef PKS_ENERGY_BASE_HH_
#define PKS_ENERGY_BASE_HH_

#include "PK_Factory.hh"

#include "PDE_Diffusion.hh"
#include "PDE_DiffusionMFD.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"

//#include "PK_PhysicalBDF_ATS.hh"
#include "pk_physical_bdf_default.hh"
#include "upwinding.hh"

namespace Amanzi {

// forward declarations
namespace Operators { class Advection; }
namespace Functions { class BoundaryFunction; }

namespace Energy {

class EnergyBase : public PK_PhysicalBDF_Default{

public:

  EnergyBase(Teuchos::ParameterList& FElist,
             const Teuchos::RCP<Teuchos::ParameterList>& plist,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& solution);
  
  // Virtual destructor
  virtual ~EnergyBase() {}

  // EnergyBase is a PK
  // -- Setup data
  virtual void Setup(const Teuchos::Ptr<State>& S) override;

  // -- Initialize owned (dependent) variables.
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) override;
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) override {}

  // Default implementations of BDFFnBase methods.
  // -- Compute a norm on u-du and return the result.
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<const TreeVector> du) override;

  // EnergyBase is a BDFFnBase
  // computes the non-linear functional f = f(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                   Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) override;


  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // problems with temperatures -- setting a range of admissible temps
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override;

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u) override;
    
  // evaluating consistent faces for given BCs and cell values
  virtual void CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u);

  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override;

 protected:
  // These must be provided by the deriving PK.
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) = 0;

  // -- get enthalpy as a function of Dirichlet boundary data.  Note that this
  //    will get replaced by a better system when we get maps on the boundary
  //    faces.
  virtual void ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S);
  virtual void ApplyDirichletBCsToBoundaryFace_(const Teuchos::Ptr<CompositeVector>& temp);

  // -- Add any source terms into the residual.
  virtual void AddSources_(const Teuchos::Ptr<State>& S,
                           const Teuchos::Ptr<CompositeVector>& f);
  virtual void AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h);

  // Standard methods
  virtual void SetupEnergy_(const Teuchos::Ptr<State>& S);

  // Upwinding conductivities
  virtual bool UpdateConductivityData_(const Teuchos::Ptr<State>& S);
  virtual bool UpdateConductivityDerivativeData_(const Teuchos::Ptr<State>& S);


  // boundary condition members
  virtual void UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S);

  // physical methods
  // -- accumulation of energy
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& f);

  // -- advection of enthalpy
  virtual void AddAdvection_(const Teuchos::Ptr<State>& S,
                     const Teuchos::Ptr<CompositeVector>& f, bool negate);

  // -- diffusion of temperature
  virtual void ApplyDiffusion_(const Teuchos::Ptr<State>& S,
          const Teuchos::Ptr<CompositeVector>& f);

  virtual int BoundaryFaceGetCell(int f) const;


 protected:
  int niter_;

  // boundary conditions
  Teuchos::RCP<Functions::BoundaryFunction> bc_temperature_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_diff_flux_;
  Teuchos::RCP<Functions::BoundaryFunction> bc_flux_;

  Teuchos::RCP<Operators::BCs> bc_adv_;

  // operators
  Teuchos::RCP<Operators::Upwinding> upwinding_;
  Teuchos::RCP<Operators::Upwinding> upwinding_deriv_;

  // mathematical operators
  Teuchos::RCP<Operators::Operator> matrix_; // pc in PKPhysicalBDFBase
  Teuchos::RCP<Operators::PDE_Diffusion> matrix_diff_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> matrix_adv_;

  Teuchos::RCP<Operators::PDE_Diffusion> preconditioner_diff_;
  Teuchos::RCP<Operators::PDE_Accumulation> preconditioner_acc_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> preconditioner_adv_;

  // flags and control
  bool modify_predictor_with_consistent_faces_;
  bool modify_predictor_for_freezing_;
  bool modify_correction_for_freezing_;
  bool is_source_term_;
  bool is_source_term_differentiable_;
  bool is_source_term_finite_differentiable_;
  bool is_mass_source_term_;
  bool implicit_advection_;
  bool implicit_advection_in_pc_;
  bool precon_used_;
  bool flux_exists_;
  bool jacobian_;

  double T_limit_;
  double mass_atol_;
  double soil_atol_;
  
  bool coupled_to_subsurface_via_temp_;
  bool coupled_to_subsurface_via_flux_;
  bool coupled_to_surface_via_temp_;
  bool coupled_to_surface_via_flux_;

  // Keys
  Key energy_key_;
  Key wc_key_;
  Key enthalpy_key_;
  Key flux_key_;
  Key energy_flux_key_;
  Key adv_energy_flux_key_;
  Key conductivity_key_;
  Key uw_conductivity_key_;
  Key dconductivity_key_;
  Key duw_conductivity_key_;
  Key source_key_;
  //  Key mass_source_key_;
  Key ss_flux_key_;
};

} // namespace Energy
} // namespace Amanzi

#endif
