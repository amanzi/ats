/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: 
*/
//! 

/*!

Solves Preferential Flow equation:


*/

#ifndef PK_FLOW_PREFERENTIAL_HH_
#define PK_FLOW_PREFERENTIAL_HH_

#include "richards.hh"

#include "PDE_DiffusionFactory.hh"
#include "PDE_Accumulation.hh"
#include "PK_Factory.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {

// forward declarations
class MPCSubsurface;
class PredictorDelegateBCFlux;
class PrimaryVariableFieldEvaluator;
namespace WhetStone { class Tensor; }

namespace Flow {

class Preferential : public Richards {

public:

  Preferential(Teuchos::ParameterList& pk_tree,
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~Preferential() {}

  // virtual void calculate_diagnostics(const Teuchos::RCP<State>& S) {CalculateDiagnostics(S);};

  // -- Setup data.
  virtual void Setup(const Teuchos::Ptr<State>& S) override;

  // // -- Initialize owned (dependent) variables.
  // virtual void Initialize(const Teuchos::Ptr<State>& S);

  // // -- Commit any secondary (dependent) variables.
  // virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);

  // // -- limit changes in a valid time step
  // virtual bool ValidStep();

  // // -- Update diagnostics for vis.
  // virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S);

  // // ConstantTemperature is a BDFFnBase
  // // computes the non-linear functional g = g(t,u,udot)
  // virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
  //                  Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // // applies preconditioner to u and returns the result in Pu
  // virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // // updates the preconditioner
  // virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
  //         Teuchos::RCP<TreeVector> u);

  // // problems with pressures -- setting a range of admissible pressures
  // virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up);

  // // evaluating consistent faces for given BCs and cell values
  // virtual void CalculateConsistentFaces(const Teuchos::Ptr<CompositeVector>& u);

protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) override;
  virtual void SetupPreferentialFlow_(const Teuchos::Ptr<State>& S);

  // // boundary condition members
  // void ComputeBoundaryConditions_(const Teuchos::Ptr<State>& S);
  // virtual void UpdateBoundaryConditions_(const Teuchos::Ptr<State>& S, bool kr=true);

  // // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  // virtual void SetAbsolutePermeabilityTensor_(const Teuchos::Ptr<State>& S);
  virtual bool UpdatePermeabilityData_(const Teuchos::Ptr<State>& S) override;
  virtual bool UpdatePermeabilityDerivativeData_(const Teuchos::Ptr<State>& S) override;

  // virtual void UpdateVelocity_(const Teuchos::Ptr<State>& S);
  // virtual void InitializeHydrostatic_(const Teuchos::Ptr<State>& S);

  // // physical methods
  // // -- diffusion term
  // virtual void ApplyDiffusion_(const Teuchos::Ptr<State>& S,
  //         const Teuchos::Ptr<CompositeVector>& g);

  // // virtual void AddVaporDiffusionResidual_(const Teuchos::Ptr<State>& S,
  // //         const Teuchos::Ptr<CompositeVector>& g);
  // // virtual void ComputeVaporDiffusionCoef(const Teuchos::Ptr<State>& S,
  // //                                        Teuchos::RCP<CompositeVector>& vapor_diff,
  // //                                        std::string var_name);

  // // -- accumulation term
  // virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);

  // // -- Add any source terms into the residual.
  // virtual void AddSources_(const Teuchos::Ptr<State>& S,
  //                          const Teuchos::Ptr<CompositeVector>& f);
  // virtual void AddSourcesToPrecon_(const Teuchos::Ptr<State>& S, double h);

  // // Nonlinear version of CalculateConsistentFaces()
  // // virtual void CalculateConsistentFacesForInfiltration_(
  // //     const Teuchos::Ptr<CompositeVector>& u);
  // virtual bool ModifyPredictorConsistentFaces_(double h, Teuchos::RCP<TreeVector> u);
  // virtual bool ModifyPredictorWC_(double h, Teuchos::RCP<TreeVector> u);
  // virtual bool ModifyPredictorFluxBCs_(double h, Teuchos::RCP<TreeVector> u);

  // // virtual void PreconWC_(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // // -- Possibly modify the correction before it is applied
  // virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  //     ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
  //                      Teuchos::RCP<const TreeVector> u,
  //                      Teuchos::RCP<TreeVector> du);

  // void  ClipHydrostaticPressure(double pmin, Epetra_MultiVector& p);

protected:
   Key coef_grav_key_, dcoef_grav_key_; 
  
private:
  // factory registration
  static RegisteredPKFactory<Preferential> reg_;

  // Preferential has a friend in couplers...
  friend class Amanzi::MPCSubsurface;

};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
