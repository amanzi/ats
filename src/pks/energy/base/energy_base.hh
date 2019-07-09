/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Base energy PK.

This provides the base of an advection-diffusion equation for energy.

------------------------------------------------------------------------- */

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
  Teuchos::RCP<Operators::Operator> lin_solver_;

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
  bool bc_surf_temp_dependent_;
};

} // namespace Energy
} // namespace Amanzi

#endif
