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


Includes options from:

* ``[pk-physical-default-spec]`` PKPhysicalDefault_
# ``[pk-bdf-default-spec]`` PKBDFDefault_
* ``[pk-physical-bdf-default-spec]`` PKPhysicalBDFDefault_

Other variable names, typically not set as the default is basically always good:

* `"conserved quantity key`" ``[string]`` **DOMAIN-water_content** Typically not set, default is good. ``[mol]``

* `"mass density key`" ``[string]`` **DOMAIN-mass_density_liquid** liquid water density ``[kg m^-3]``

* `"molar density key`" ``[string]`` **DOMAIN-molar_density_liquid** liquid water density ``[mol m^-3]``

* `"permeability key`" ``[string]`` **DOMAIN-permeability** permeability of the soil medium ``[m^2]``

* `"conductivity key`" ``[string]`` **DOMAIN-relative_permeability** scalar coefficient of the permeability ``[-]``

* `"upwind conductivity key`" ``[string]`` **DOMAIN-upwind_relative_permeability** upwinded (face-based) scalar coefficient of the permeability.  Note the units of this are strange, but this represents :math:`\frac{n_l k_r}{\mu}`  ``[mol kg^-1 s^1 m^-2]``

* `"darcy flux key`" ``[string]`` **DOMAIN-mass_flux** mass flux across a face ``[mol s^-1]``

* `"darcy flux direction key`" ``[string]`` **DOMAIN-mass_flux_direction** direction of the darcy flux (used in upwinding :math:`k_r`) ``[??]``

* `"darcy velocity key`" ``[string]`` **DOMAIN-darcy_velocity** darcy velocity vector, interpolated from faces to cells ``[m s^-1]``

* `"darcy flux key`" ``[string]`` **DOMAIN-mass_flux** mass flux across a face ``[mol s^-1]``

* `"saturation key`" ``[string]`` **DOMAIN-saturation_liquid** volume fraction of the liquid phase ``[-]``

Discretization control:

* `"diffusion`" ``[list]`` An PDE_Diffusion_ spec describing the (forward) diffusion operator

* `"diffusion preconditioner`" ``[list]`` An PDE_Diffusion_ spec describing the diffusive parts of the preconditioner.

* `"linear solver`" ``[linear-solver-typed-spec]`` **optional** is a LinearSolver_ spec.  Note
  that this is only used if this PK is not strongly coupled to other PKs.

Boundary conditions:

//* `"boundary conditions`" ``[subsurface-flow-bc-spec]`` Defaults to Neuman, 0 normal flux.  See `Flow-specific Boundary Conditions`_

Physics control:

* `"permeability rescaling`" ``[double]`` **1** Typically 1e7 or order :math:`sqrt(K)` is about right.  This rescales things to stop from multiplying by small numbers (permeability) and then by large number (:math:`\rho / \mu`).

* `"permeability type`" ``[string]`` **'scalar'** The permeability type can be 'scalar', 'horizontal and vertical', 'diagonal tensor', or 'full tensor'. This key is placed in state->field evaluators->permeability. The 'scalar' option requires 1 permeability value, 'horizontal and vertical' requires 2 values, 'diagonal tensor' requires 2 (2D) or 3 (3D) values, and 'full tensor' requires 3 (2D) or 6 (3D) values. The ordering of the permeability values in the input script is important: 'horizontal and vertical'={xx/yy,zz}, 'diagonal tensor'={xx,yy} or {xx,yy,zz}, 'full tensor'={xx,yy,xy/yx} or {xx,yy,zz,xy/yx,xz/zx,yz/zy}.

* `"water retention evaluator`" ``[wrm-evaluator-spec]`` The WRM.  This needs to go away!

This PK additionally requires the following:

EVALUATORS:
- `"conserved quantity`"
- `"mass density`"
- `"molar density`"
- `"permeability`"
- `"conductivity`"
- `"saturation`"

*/

#pragma once

#include "errors.hh"

#include "TreeVector.hh"
#include "SolverDefs.hh"

#include "PK_Adaptors.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinLeaf.hh"
#include "PK_Default.hh"
#include "PK_Factory.hh"

namespace ATS {
namespace Flow {

using namespace Amanzi;

template <class Base_t>
class PK_Richards : public Base_t {

public:

  PK_Richards(const Teuchos::RCP<Teuchos::ParameterList>& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
              const Teuchos::RCP<State>& S)
    : Base_t(pk_tree, global_plist, S)
  {}

  // -- Setup data.
  void Setup() {
    Base_t::Setup();

    // require evaluator
    S_->RequireEvaluator(key_, tag_inter_);
    
    S_->template Require<CompositeVector, CompositeVectorSpace>(key_, tag_old_, key_)
        .SetMesh(mesh_)
        ->SetComponent("cell", AmanziMesh::CELL, 1);

    S_->template Require<CompositeVector, CompositeVectorSpace>(key_, tag_new_, key_)
        .SetMesh(mesh_)
        ->SetComponent("cell", AmanziMesh::CELL, 1);
        
  }

  // -- Initialize owned (primary) variables.
  void Initialize() {
    Base_t::Initialize();
    auto& ic_cv = S_->template GetW<CompositeVector>(key_, "", key_);
    ic_cv.putScalar(0.);
  }

  void FunctionalTimeDerivative(double t, const TreeVector& u, TreeVector& f) {
    f.putScalar(0.);
  }

  void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                     Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) {
    double dt = t_new - t_old;
    db_->WriteVector("u_old", *u_old->Data());
    db_->WriteVector("u_new", *u_new->Data());

    f->update(1/dt, *u_new, -1/dt, *u_old, 0.);
    db_->WriteVector("res", *f->Data());
  }

  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                          Teuchos::RCP<TreeVector> Pu) {
    *Pu = *u;
    return 0;
  }

  void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {}

  double
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) {
    return du->normInf();
  }

  void ChangedSolution() { this->ChangedSolutionPK(tag_new_); }

  void UpdateContinuationParameter(double lambda) {
    Errors::Message msg("Richards PK: continuation methods not supported.");
    throw(msg);
  }
  
  // -- limit changes in a valid time step
  bool ValidStep(const Key& tag_old, const Key& tag_new) {
    return true;
  }

  // -- Update diagnostics for vis.
  void CalculateDiagnostics(const Key& tag) {}

  bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                       Teuchos::RCP<TreeVector> u) {
    return false;
  }

  // problems with pressures -- setting a range of admissible pressures
  bool IsAdmissible(Teuchos::RCP<const TreeVector> up) {
    return true;
  }

  // -- Possibly modify the correction before it is applied
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du){
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

 protected:
  using Base_t::tag_new_;
  using Base_t::tag_old_;
  using Base_t::tag_inter_;
  using Base_t::S_;
  using Base_t::key_;
  using Base_t::mesh_;
  using Base_t::db_;

  Key dudt_key_;
  Key res_key_;
  Key conserved_key_;
  
 private:
};

using Richards_Implicit =
    PK_Implicit_Adaptor<PK_Richards<
                          PK_MixinImplicit<
                            PK_MixinLeafCompositeVector<PK_Default>>>>;

using Richards_Explicit =
    PK_Explicit_Adaptor<PK_Richards<
                          PK_MixinExplicit<
                            PK_MixinLeafCompositeVector<PK_Default>>>>;

    
    
}  // namespace Flow
}  // namespace ATS

