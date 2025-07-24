/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/
/*!

This MPC provides terms for an approximate Jacobian for coupling the two-phase
diffusion wave equation (the :ref:`Overland Flow with Ice`) to the two-phase
surface energy equation (the :ref:`Overland Energy with Ice`).

`"PK type`" = `"icy surface`"

.. _pk-icy-surface-spec:
.. admonition:: pk-icy-surface-spec

   * `"domain name`" ``[string]`` Domain of simulation

   * `"preconditioner type`" ``[string]`` **picard** See the above for detailed
     descriptions of the choices.  One of: `"none`", `"block diagonal`", `"no
     flow coupling`", `"picard`", `"ewc`", and `"smart ewc`".  Note that the
     `"ewc`" preconditioners are rarely helpful in integrated problems, but may
     be helpful in surface-only problems.

   * `"supress Jacobian terms: d div surface q / dT`" ``[bool]`` **false** Turn
     off the Jacobian term associated with the derivative of overland
     conductivity with respect to temperature.

   * `"surface ewc delegate`" ``[mpc-delegate-ewc-spec]`` A :ref:`EWC Globalization
     Delegate` spec.

   KEYS

   - `"temperature`"
   - `"pressure`"
   - `"energy`"
   - `"water content`"
   - `"overland conductivity`"
   - `"upwind overland conductivity`"
   - `"potential`"
   - `"ponded depth, negative`"
   - `"water flux`"

*/

#ifndef MPC_SURFACE_HH_
#define MPC_SURFACE_HH_

#include "TreeOperator.hh"
#include "pk_physical_bdf_default.hh"
#include "strong_mpc.hh"

namespace Amanzi {
namespace Operators {
class PDE_Diffusion;
class PDE_Accumulation;
class Operator;
class UpwindTotalFlux;
class UpwindArithmeticMean;
class Upwinding;
} // namespace Operators
class MPCDelegateEWCSurface;

class MPCSurface : public StrongMPC<PK_PhysicalBDF_Default> {
 public:
  MPCSurface(Teuchos::ParameterList& pk_tree_list,
             const Teuchos::RCP<Teuchos::ParameterList>& global_list,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& soln);

  void parseParameterList() override;
  virtual void Setup() override;
  virtual void Initialize() override;
  //virtual void set_tags(const Tag& tag_current, const Tag& tag_next);

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h,
                               Teuchos::RCP<const TreeVector> up0,
                               Teuchos::RCP<TreeVector> up) override;

  // -- Update the preconditioner to be physically consistent
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> Pu) override;

  Teuchos::RCP<Operators::TreeOperator> preconditioner() { return preconditioner_; }

 protected:
  enum PreconditionerType {
    PRECON_NONE = 0,
    PRECON_BLOCK_DIAGONAL = 1,
    PRECON_PICARD = 2,
    PRECON_EWC = 3,
    PRECON_NO_FLOW_COUPLING = 4,
  };

  Teuchos::RCP<Operators::TreeOperator> preconditioner_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  // preconditioner methods
  PreconditionerType precon_type_;

  // Additional precon terms
  //   equations are given by:
  // 1. conservation of WC: dWC/dt + div q = 0
  // 2. conservation of E:  dE/dt + div K grad T + div hq = 0

  // dWC / dT off-diagonal block
  Teuchos::RCP<Operators::Operator> dWC_dT_block_;
  // -- d ( div q ) / dT  terms
  Teuchos::RCP<Operators::PDE_Diffusion> ddivq_dT_;

  // dE / dp off-diagonal block
  Teuchos::RCP<Operators::Operator> dE_dp_block_;
  // -- d ( div hq ) / dp terms
  //Teuchos::RCP<Operators::PDE_Diffusion> ddivhq_dp_;
  // -- d ( dE/dt ) / dp terms
  Teuchos::RCP<Operators::PDE_Accumulation> dE_dp_;

  // dE / dT on-diagonal block additional terms that use q info
  // -- d ( div hq ) / dT terms
  //Teuchos::RCP<Operators::PDE_Diffusion> ddivhq_dT_;

  Key domain_;
  Key temp_key_;
  Key pres_key_;
  Key e_key_;
  Key wc_key_;
  Key tc_key_;
  Key kr_key_;
  Key kr_uw_key_;
  Key potential_key_;
  Key pd_bar_key_;
  Key enth_key_;
  Key water_flux_key_;
  Key rho_key_;

  // EWC delegate
  Teuchos::RCP<MPCDelegateEWCSurface> ewc_;

  // cruft for easier global debugging
  bool dump_;
  Teuchos::RCP<Debugger> db_;

  // do we need to rescale preconditioner to pressure from head?
  bool rescale_precon_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSurface> reg_;
};


} // namespace Amanzi

#endif
