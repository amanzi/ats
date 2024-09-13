/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

//! A coupler which integrates surface, richards and preferantial flows implicitly.
/*!



*/

#ifndef PKS_MPC_COUPLED_DUALMEDIA_WATER_HH_
#define PKS_MPC_COUPLED_DUALMEDIA_WATER_HH_

#include "Operator.hh"
#include "mpc_delegate_water.hh"
#include "pk_physical_bdf_default.hh"
#include "TreeOperator.hh"
#include "strong_mpc.hh"

namespace Amanzi {

class MPCCoupledDualMediaWater : public StrongMPC<PK_BDF_Default> {
 public:
  MPCCoupledDualMediaWater(Teuchos::ParameterList& FElist,
                           const Teuchos::RCP<Teuchos::ParameterList>& plist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& soln);

  void parseParameterList() override;

  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk FunctionalResidual().
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g);

  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  // -- Update the preconditioner.
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

 protected:
  // void
  // UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
  //         const Teuchos::RCP<TreeVector>& Pu);

 protected:
  Teuchos::RCP<Operators::TreeOperator> op_tree_matrix_, op_tree_pc_;
  Teuchos::RCP<TreeVector> op_tree_rhs_;

  // sub PKs
  Teuchos::RCP<PK_BDF_Default> surf_flow_pk_;
  Teuchos::RCP<PK_BDF_Default> macro_flow_pk_;
  Teuchos::RCP<StrongMPC<PK_PhysicalBDF_Default>> integrated_flow_pk_;
  Teuchos::RCP<PK_BDF_Default> matrix_flow_pk_;


  Key domain_;
  Key domain_macropore_;
  Key total_ss_flux_key_, matrix_flux_key_, macro_flux_key_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> macropore_db_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCCoupledDualMediaWater> reg_;
};

} // namespace Amanzi


#endif
