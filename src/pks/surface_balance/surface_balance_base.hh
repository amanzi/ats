/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A simple conservation ODE.

/*!

This is a very simple vector of ODEs, useful in balance equations, where the
time derivative of a conserved quantity is determined by a bunch of sources and
sinks.

.. math::
    \frac{\partial \Phi }{\partial t} = \sum_i Q_i

.. _balance-pk-spec:
.. admonition:: balance-pk-spec

    * `"domain`" ``[string]`` Mesh on which the balance is to be done.

    * `"primary variable key`" ``[string]`` The primary variable associated with
      this PK.  Note there is no default -- this must be provided by the user.

    * `"conserved quantity key`" ``[string]`` The conserved quantity :math:`\Phi`

    * `"source key`" ``[string]`` **DOMAIN-source_sink** Units are in conserved
      quantity per second per cell volume.

    * `"time discretization theta`" ``[double]`` **1.0** :math:`\theta` in a
      Crank-Nicholson time integration scheme.  1.0 implies fully implicit, 0.0
      implies explicit, 0.5 implies C-N.

    * `"modify predictor positivity preserving`" ``[bool]`` **false** If true,
      predictors are modified to ensure that the conserved quantity is always > 0.

    * `"absolute error tolerance`" ``[double]`` **550.0** a_tol in the standard
      error norm calculation.  Defaults to a small amount of water.  Units are
      the same as the conserved quantity.

    * `"preconditioner`" ``[preconditioner-typed-spec]`` Preconditioner for the
      solve.

    * `"linear solver`" ``[linear-solver-typed-spec]`` **optional** May be used
      to improve the inverse of the diffusion preconditioner.  Only used if this
      PK is not implicitly coupled.  See LinearOperator_.

    INCLUDES:

    - ``[pk-physical-bdf-default-spec]``
      
*/

#ifndef PK_SURFACE_BALANCE_BASE_HH_
#define PK_SURFACE_BALANCE_BASE_HH_

#include "Operator.hh"
#include "PDE_Accumulation.hh"
#include "PK_Factory.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {
namespace SurfaceBalance {

class SurfaceBalanceBase : public PK_PhysicalBDF_Default {

 public:

  SurfaceBalanceBase(Teuchos::ParameterList& pk_tree,
                     const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                     const Teuchos::RCP<State>& S,
                     const Teuchos::RCP<TreeVector>& solution);

  // main methods
  // -- Setup data.
  virtual void Setup(const Teuchos::Ptr<State>& S) override;

  // -- Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) override {}

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> g) override;

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override;

 protected:

  Key layer_;
  bool conserved_quantity_;
  bool is_source_, is_source_differentiable_, source_finite_difference_;
  Key source_key_;

  double theta_;
  double eps_;

  bool modify_predictor_positivity_preserving_;

  bool precon_used_;
  Teuchos::RCP<Operators::PDE_Accumulation> preconditioner_acc_;
  Teuchos::RCP<Operators::Operator> lin_solver_;
  
 private:
  // factory registration
  static RegisteredPKFactory<SurfaceBalanceBase> reg_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif

