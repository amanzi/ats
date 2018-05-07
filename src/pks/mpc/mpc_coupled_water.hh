/* -*-  mode: c++; indent-tabs-mode: nil -*- */
//! MPCCoupledWater: coupler which integrates surface and subsurface flow.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


#ifndef PKS_MPC_COUPLED_WATER_HH_
#define PKS_MPC_COUPLED_WATER_HH_

/*!

Couples Richards equation to surface water through continuity of both pressure and fluxes.

Currently requires that the subsurface discretization is a face-based
discretization, i.e. one of the MFD methods.  Then the surface equations are
directly added into the subsurface discrete equations.

* `"PKs order`" ``[Array(string)]`` Supplies the names of the coupled PKs.
  The order must be {subsurface_flow_pk, surface_flow_pk} (subsurface first).

* `"linear solver`" ``[linear-solver-spec]`` A LinearSolver_ spec.  Only used
  if this PK is not itself coupled by other strong couplers.

* `"preconditioner`" ``[preconditioner-spec]`` A Preconditioner_ spec.  Only used
  if this PK is not itself coupled by other strong couplers.

* `"water delegate`" ``[list]`` 

*/


#include "Operator.hh"
#include "mpc_delegate_water.hh"
#include "pk_physical_bdf_default.hh"
#include "LinearOperatorFactory.hh"
#include "strong_mpc.hh"

namespace Amanzi {

class MPCCoupledWater : public StrongMPC<PK_PhysicalBDF_Default> {
 public:


  MPCCoupledWater(Teuchos::ParameterList& FElist,
                  const Teuchos::RCP<Teuchos::ParameterList>& plist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk Functional().
  virtual void Functional(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
           Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g);

  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
          Teuchos::RCP<TreeVector> u);

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> du);

 protected:
  // void
  // UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
  //         const Teuchos::RCP<TreeVector>& Pu);

 protected:

  // sub PKs
  Teuchos::RCP<PK_PhysicalBDF_Default> domain_flow_pk_;
  Teuchos::RCP<PK_PhysicalBDF_Default> surf_flow_pk_;

  // sub meshes
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;

  // coupled preconditioner
  Teuchos::RCP<Operators::Operator> precon_;
  Teuchos::RCP<Operators::Operator> precon_surf_;
  Teuchos::RCP<Operators::Operator> lin_solver_;

  // Water delegate
  Teuchos::RCP<MPCDelegateWater> water_;
  bool consistent_cells_;
  
  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCCoupledWater> reg_;

  std::string domain_surf, domain_ss;

};

} // namespace


#endif
