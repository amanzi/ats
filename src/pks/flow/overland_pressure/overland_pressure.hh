/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
This is the overland flow component of ATS.
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
----------------------------------------------------------------------------- */

#ifndef PK_FLOW_OVERLAND_HEAD_HH_
#define PK_FLOW_OVERLAND_HEAD_HH_

#include "BoundaryFunction.hh"
#include "upwinding.hh"

#include "Operator.hh"
#include "OperatorDiffusion.hh"
#include "OperatorAccumulation.hh"

#include "pk_factory.hh"
#include "pk_physical_bdf_base.hh"

namespace Amanzi {

class MPCSurfaceSubsurfaceDirichletCoupler;
namespace Flow {

namespace FlowRelations {
  class OverlandConductivityModel;
  class HeightModel;
}


class OverlandPressureFlow : public PKPhysicalBDFBase {

public:
  OverlandPressureFlow(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<TreeVector>& solution);
  
  // Virtual destructor
  virtual ~OverlandPressureFlow() {}

  // main methods
  // -- Initialize owned (dependent) variables.
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- Initialize owned (dependent) variables.
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // -- Commit any secondary (dependent) variables.
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);

  // -- Update diagnostics for vis.
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
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

  virtual void FixBCsForOperator_(const Teuchos::Ptr<State>& S);
  virtual void FixBCsForPrecon_(const Teuchos::Ptr<State>& S);
  // virtual void FixBCsForConsistentFaces_(const Teuchos::Ptr<State>& S);

  virtual void ApplyBoundaryConditions_(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<CompositeVector>& pres );

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
  bool source_in_meters_;
  bool modify_predictor_with_consistent_faces_;
  bool symmetric_;
  bool perm_update_required_;
  bool source_only_if_unfrozen_;
  bool smoothed_ponded_accumulation_;

  double p_limit_;
  double patm_limit_;

  // coupling term
  bool coupled_to_subsurface_via_head_;
  bool coupled_to_subsurface_via_flux_;
  Key mass_source_key_;
  bool jacobian_;

  // work data space
  Teuchos::RCP<Operators::Upwinding> upwinding_;
  Teuchos::RCP<Operators::Upwinding> upwinding_dkdp_;

  // mathematical operators
  Teuchos::RCP<Operators::Operator> matrix_; // pc in PKPhysicalBDFBase
  Teuchos::RCP<Operators::OperatorDiffusion> matrix_diff_;
  Teuchos::RCP<Operators::OperatorDiffusion> face_matrix_diff_;
  Teuchos::RCP<Operators::OperatorDiffusion> preconditioner_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> preconditioner_acc_;
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
  
  // needed physical models
  Teuchos::RCP<FlowRelations::OverlandConductivityModel> cond_model_;

  int niter_;

  // factory registration
  static RegisteredPKFactory<OverlandPressureFlow> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
