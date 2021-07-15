/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A parabolic PDE in mixed form for conservation.

/*!

This is the canonical nonlinear parabolic PDE, in mixed form.

.. math::
    \frac{\partial \Psi(u) }{\partial t} - \nabla \cdot K(u) \grad \Phi(u) = Q(u,x,t)

Where the conserved quantity :math:`\Psi` is a function of the primary variable
:math:`u`, diffusive fluxes are provided as a function of the coefficient
:math:`K` and gradients in a potential :math:`\Phi` which is also a function of
the primary variable, and a source term :math:`Q` is provided.

Note that while this is mixed form, some classes here may solve it in the
primary form,

.. math::
    \frac{d \Psi(u) }{d u} \frac{\partial u}{\partial t} - \nabla \cdot K(u) \grad \Phi(u) = Q(u,x,t)

where we then must assume that :math:`\frac{d \Psi(u) }{d u} > 0`.

Note that frequently \Phi(u) = u (e.g. the potential field is the primary
variable) but not always -- specifically for surface water where the potential
field is :math:`h(p) + z` for primary variable pressure.

.. _conservation-ode-pk-spec:
.. admonition:: conservation-ode-pk-spec

    * `"domain`" ``[string]`` Mesh on which the balance is to be done.

    * `"primary variable key`" ``[string]`` The primary variable, :math:`u`.
      Note there is no default -- this must be provided by the user.

    * `"conserved quantity key`" ``[string]`` The conserved quantity :math:`\Psi`

    * `"diffusion operand key`" ``[string]`` The diffused quantity :math:`\Phi`

    * `"source key`" ``[string]`` **DOMAIN-source_sink** Units are in conserved
      quantity per second, :math:`Q`.

    * `"time discretization theta`" ``[double]`` **1.0** :math:`\theta` in a
      Crank-Nicholson time integration scheme.  1.0 implies fully implicit, 0.0
      implies explicit, 0.5 implies C-N.  Note, only used in the implicit
      scheme -- prefer to use an explicit time integrator over :math:`\theta ==
      0`.

    * `"modify predictor positivity preserving`" ``[bool]`` **false** If true,
      predictors are modified to ensure that the conserved quantity is always >
      0.  These systems may be stiff, and this does not guarantee positivity,
      so time integration methods may need to be chosen with care.

    * `"absolute error tolerance`" ``[double]`` **550.0** a_tol in the standard
      error norm calculation.  Defaults to a small amount of water.  Units are
      the same as the conserved quantity (default assumes mols).

    * `"inverse`" ``[inverse-typed-spec]`` **optional**
      Linear inverse for the linear solve.  Only used if the time integration
      scheme is solved implicitly.

*/

#pragma once

#include "TreeVector.hh"
#include "Operator.hh"
#include "Operator_Factory.hh"
#include "Inverse.hh"
#include "SolverDefs.hh"

#include "PK_Adaptors.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinPredictorCorrector.hh"
#include "PK_MixinLeaf.hh"
#include "PK_Default.hh"
#include "PK_Factory.hh"

#include "PK_MixinConservationEquation.hh"

namespace ATS {
namespace Basic {

using namespace Amanzi;

template <class Base_t>
class MixedFormParabolicPDE_Implicit : public Base_t {

  public:
  using Base_t::Base_t;

  void Setup()
  {
    // Call any other mixin Setup methods, including time integration setup.
    Base_t::Setup();

    // Some mixins need a time tag at which to require their evaluators
    // (e.g. conserved quantity, etc).
    Base_t::SetupAtTag(tag_new_);

    // Set up accumulation
    // -------------------
    // This evaluator calculates dPsi/dt, a function of Psi and time at two
    // time tags.
    // -- Note the suffix _t indicates a time derivative -- the output
    accumulation_key_ = conserved_quantity_key_ + "_t";
    Teuchos::ParameterList& acc_list = S_->FEList().sublist(accumulation_key_);
    // -- Note, I'm not sure why a user would override these by setting them in
    //    the parameterlist manually, but have chosen not to error or overwrite
    //    the user's values in case that changes.
    if (!acc_list.isParameter("conserved quantity key"))
      acc_list.set("conserved quantity key", conserved_quantity_key_);
    // -- Note that by not hard-coding tag_old and tag_new, this same PK is used
    //    for subcycling.
    if (!acc_list.isParameter("tag")) acc_list.set("tag", tag_new_);
    if (!acc_list.isParameter("tag old")) acc_list.set("tag old", tag_old_);
    if (!acc_list.isParameter("tag new")) acc_list.set("tag new", tag_new_);
    if (!acc_list.isParameter("evaluator type")) acc_list.set("evaluator type", "accumulation operator");

    // -- Require a vector to store the output dPsi/dt.
    S_->template Require<CompositeVector,CompositeVectorSpace>(accumulation_key_, tag_new_)
        .SetMesh(mesh_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);

    // Set up diffusion
    // ----------------
    // This evaluator calculates q and local matrices for MFD or FV
    // representations of the div k K grad operator.  The evaluator itself is
    // the PDE_Diffusion(withGravity) object.
    //
    // Note we don't actually use the local matrices in this PK, so we don't
    // require them.  We require the residual evaluator (which is a global
    // operator) and that in turn requires the local matrices.  But by setting
    // up the parameter list here, we can make the input file simpler and save
    // the user from having to write all this themselves.
    //
    // -- Input includes a diffused variable -- allows primary and diffused variables to be different
    diffused_key_ = Keys::readKey(*plist_, domain_, "diffusion operand", Keys::getVarName(key_));
    if (diffused_key_ != key_)
      S_->template RequireDerivative<CompositeVector,CompositeVectorSpace>(diffused_key_, tag_new_, key_, tag_new_);
    // -- Require the diffusion evaluator (a PDE_Diffusion object)
    diffusion_key_ = Keys::readKey(*plist_, domain_, "diffusion operator");
    Teuchos::ParameterList& diff_list = S_->FEList().sublist(diffusion_key_);
    if (!diff_list.isParameter("evaluator type")) diff_list.set("evaluator type", "diffusion operator");
    if (!diff_list.isParameter("local operator key")) diff_list.set("local operator key", diffusion_key_);
    if (!diff_list.isParameter("rhs key")) diff_list.set("rhs key", diffusion_key_+"_rhs");
    if (!diff_list.isParameter("flux key")) diff_list.set("flux key", diffusion_key_+"_flux");
    if (!diff_list.isParameter("boundary conditions key")) diff_list.set("boundary conditions key", diffusion_key_+"_bcs");
    diff_list.set("operator argument key", diffused_key_);

    // Set up the residual evaluator
    // -----------------------------
    // This evaluator is the global operator (or potentially a collection of
    // global operators), and calculates r = A*p - b
    res_key_ = key_ + "_res";
    Teuchos::ParameterList& res_list = S_->FEList().sublist(res_key_);
    // -- operator at the new time
    if (!res_list.isParameter("tag")) res_list.set("tag", tag_new_);
    // -- type is an Evaluator_OperatorApply, which uses an Operator object.
    if (!res_list.isParameter("evaluator type")) res_list.set("evaluator type", "operator application");
    // -- the global operand on the diagonal (for coupled PKs, the MPC could
    //    add additional terms for offdiagonal operators)
    res_list.set("diagonal primary x key", diffused_key_);
    if (!res_list.isParameter("diagonal local operators keys")) res_list.set("diagonal local operators keys",
                 Teuchos::Array<std::string>(1, diffusion_key_));
    if (!res_list.isParameter("diagonal local operator rhss keys")) res_list.set("diagonal local operator rhss keys",
                 Teuchos::Array<std::string>(1, diffusion_key_+"_rhs"));

    // -- rhs for source and accumulation term
    Teuchos::Array<std::string> rhss;
    Teuchos::Array<double> rhs_coefs;
    if (is_source_) {
      rhss = std::vector<std::string>{ source_key_, accumulation_key_ };
      rhs_coefs = std::vector<double>{ -1.0, 1.0 };
    } else {
      rhss = std::vector<std::string>{ accumulation_key_ };
      rhs_coefs = std::vector<double>{ 1.0 };
    }
    if (!res_list.isParameter("additional rhss keys")) res_list.set("additional rhss keys", rhss);
    if (!res_list.isParameter("rhs coefficients")) res_list.set("rhs coefficients", rhs_coefs);

    // -- preconditioner -- could be provided in this PK's list, or in the
    //    State evaluator list.
    if (!res_list.isSublist("inverse")) res_list.set("inverse", plist_->sublist("inverse"));


    // NOTE: we cannot know the structure of u here -- it may be CELL if FV, or
    // it may be CELL+FACE if MFD.  It will get set by the operator.  But we do
    // need to supply the mesh.
    S_->template Require<CompositeVector,CompositeVectorSpace>(key_, tag_new_)
        .SetMesh(mesh_);
    S_->template Require<CompositeVector,CompositeVectorSpace>(res_key_, tag_new_)
        .SetMesh(mesh_);
    S_->template RequireDerivative<Operators::Operator, Operators::Operator_Factory>(
      res_key_, tag_new_, diffused_key_, tag_new_);
    S_->RequireEvaluator(res_key_, tag_new_);
  }


  void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                     Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f)
  {
    // implicit time integration boilerplate
    // NOTE: Currently, we do not allow time integrators to introduce their own
    // tags -- time integrators have no knowledge of tags.  This makes it easy
    // for time integrators to break our state model, by simply
    // copy-constructing the input vector and giving us the copy instead of the
    // vector at the right tag.
    //
    // Implicit time integrators are a bit better than explicit, as they know
    // about ChangedSolution(), but better safe than sorry.  For now, these
    // check to ensure that the time integrator is giving us what we expect,
    // and is playing nice with state
    AMANZI_ASSERT(S_->template GetPtr<CompositeVector>(key_, tag_old_).get() == u_old->Data().get());
    AMANZI_ASSERT(S_->template GetPtr<CompositeVector>(key_, tag_new_).get() == u_new->Data().get());
    AMANZI_ASSERT(std::abs(S_->time(tag_old_) - t_old) < 1.e-6);
    AMANZI_ASSERT(std::abs(S_->time(tag_new_) - t_new) < 1.e-6);
    // end boilerplate

    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "----------------------------------------------------------------" << std::endl
                 << "Residual calculation: t0 = " << t_old
                 << " t1 = " << t_new << " dt = " << t_new - t_old << std::endl;

    db_->WriteCellInfo(true);
    db_->WriteVectors({key_+"_old",key_+"_new"},
                      {u_old->Data().ptr(), u_new->Data().ptr()},
                      true);

    S_->GetEvaluator(res_key_, tag_new_).Update(*S_, this->name());
    f->Data()->assign(S_->template Get<CompositeVector>(res_key_, tag_new_));

    db_->WriteState(*S_, tag_old_);
    db_->WriteState(*S_, tag_new_);
    db_->WriteVector(res_key_, f->Data().ptr(), true);
  }

  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> r,
                          Teuchos::RCP<TreeVector> Pr)
  {
    const auto& lin_op = S_->template GetDerivativePtr<Operators::Operator>(res_key_, tag_new_, diffused_key_, tag_new_);

    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "Precon application:" << std::endl;

    db_->WriteVector("u_res", r->Data().ptr());
    int ierr = lin_op->applyInverse(*r->Data(), *Pr->Data());
    db_->WriteVector("PC*u_res", Pr->Data().ptr());

    if (diffused_key_ != key_) {
      // primary and diffused key are NOT the same, must now apply the chain rule
      const auto& ddiffused_dprimary = S_->template GetDerivative<CompositeVector>(diffused_key_, tag_new_, key_, tag_new_);
      // *vo_->os() << "d " << diffused_key_ << " d " << key_ << std::endl;
      // ddiffused_dprimary.Print(*vo_->os());
      // NOTE: this is sloppy and likely needs to be put in a kernel
      Pr->Data()->reciprocal(*Pr->Data());
      Pr->Data()->elementWiseMultiply(1.0, ddiffused_dprimary, *Pr->Data(), 0.);
      Pr->Data()->reciprocal(*Pr->Data());
      // db_->WriteVector("dprimary/ddiffused*PC*u_res", Pr->Data().ptr());
    }

    // return ierr;
    return 0;  // keep trying, even if linear solver fails to converge.
  }

  void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "Precon update at t = " << t << std::endl;

    S_->GetEvaluator(res_key_, tag_new_).UpdateDerivative(*S_, this->name(), diffused_key_, tag_new_);
    if (diffused_key_ != key_)
      S_->GetEvaluator(diffused_key_, tag_new_).UpdateDerivative(*S_, this->name(), key_, tag_new_);
  }

  // -- Update diagnostics for vis.
  void CalculateDiagnostics(const Key& tag) {}

 protected:
  using Base_t::plist_;
  using Base_t::domain_;
  using Base_t::tag_new_;
  using Base_t::tag_old_;
  using Base_t::S_;
  using Base_t::key_;
  using Base_t::mesh_;
  using Base_t::db_;
  using Base_t::vo_;
  using Base_t::conserved_quantity_key_;
  using Base_t::is_source_;
  using Base_t::source_key_;

  Key accumulation_key_;
  Key diffused_key_;
  Key diffusion_key_;
  Key res_key_;

};


template <class Base_t>
class MixedFormParabolicPDE_Explicit : public Base_t {

public:
  using Base_t::Base_t;

  void Setup()
  {
    Base_t::Setup();
    Base_t::SetupAtTag(tag_inter_);

    // set up the diffusion operator
    diffusion_key_ = Keys::readKey(*plist_, domain_, "diffusion operator");
    Teuchos::ParameterList& diff_list = S_->FEList().sublist(diffusion_key_);
    if (!diff_list.isParameter("evaluator type")) diff_list.set("evaluator type", "diffusion operator");
    if (!diff_list.isParameter("local operator key")) diff_list.set("local operator key", diffusion_key_);
    if (!diff_list.isParameter("rhs key")) diff_list.set("rhs key", diffusion_key_+"_rhs");
    if (!diff_list.isParameter("boundary conditions key")) diff_list.set("boundary conditions key", diffusion_key_+"_bcs");
    if (!diff_list.isParameter("operator argument key")) diff_list.set("operator argument key", key_);

    S_->template Require<CompositeVector,CompositeVectorSpace>(conserved_quantity_key_, tag_inter_)
        .SetMesh(mesh_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
    S_->template RequireDerivative<CompositeVector,CompositeVectorSpace>(conserved_quantity_key_, tag_inter_,
            key_, tag_inter_);
    S_->RequireEvaluator(conserved_quantity_key_, tag_inter_);

    // require du/dt
    dudt_key_ = conserved_quantity_key_ + "_t";
    Teuchos::ParameterList& res_list = S_->FEList().sublist(dudt_key_);
    // -- operator
    if (!res_list.isParameter("tag")) res_list.set("tag", tag_inter_);
    if (!res_list.isParameter("evaluator type")) res_list.set("evaluator type", "operator application");
    if (!res_list.isParameter("diagonal primary x key")) res_list.set("diagonal primary x key", key_);
    if (!res_list.isParameter("diagonal local operators keys")) res_list.set("diagonal local operators keys",
                 Teuchos::Array<std::string>(1, diffusion_key_));
    if (!res_list.isParameter("diagonal local operator rhss keys")) res_list.set("diagonal local operator rhss keys",
                 Teuchos::Array<std::string>(1, diffusion_key_+"_rhs"));

    // -- rhs for source and accumulation term
    Teuchos::Array<std::string> rhss;
    Teuchos::Array<double> rhs_coefs;
    if (is_source_) {
      rhss = std::vector<std::string>{ source_key_ };
      rhs_coefs = std::vector<double>{ -1.0 };
    }
    if (!res_list.isParameter("additional rhss keys")) res_list.set("additional rhss keys", rhss);
    if (!res_list.isParameter("rhs coefficients")) res_list.set("rhs coefficients", rhs_coefs);

    // NOTE: we cannot know the structure of u here -- it may be CELL if FV, or
    // it may be CELL+FACE if MFD.  It will get set by the operator.  But we do
    // need to supply the mesh.
    S_->template Require<CompositeVector,CompositeVectorSpace>(dudt_key_, tag_inter_)
        .SetMesh(mesh_);
    S_->RequireEvaluator(dudt_key_, tag_inter_);

  }

  void
  FunctionalTimeDerivative(double t, const TreeVector& u, TreeVector& f)
  {
    // explicit time integration boilerplate
    // NOTE: Currently, we do not allow time integrators to introduce their own
    // tags -- time integrators have no knowledge of tags.  This makes it easy
    // for time integrators to break our state model, by simply
    // copy-constructing the input vector and giving us the copy instead of the
    // vector at the right tag.
    //
    // For now, these check to ensure that the time integrator is giving
    // us what we expect, and is playing nice with state
    S_->set_time(tag_inter_, t);
    AMANZI_ASSERT(u.Data() == S_->template GetPtr<CompositeVector>(key_, tag_inter_));
    this->ChangedSolutionPK(tag_inter_);
    // end boilerplate

    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "----------------------------------------------------------------" << std::endl
                 << "Time derivative calculation at: t = " << t << std::endl;

    db_->WriteCellInfo(true);
    db_->WriteVector("u", u.Data().ptr());

    S_->GetEvaluator(dudt_key_, tag_inter_).Update(*S_, this->name());
    db_->WriteVector("-div grad u - Q", S_->template GetPtr<CompositeVector>(dudt_key_, tag_inter_).ptr());

    S_->GetEvaluator(conserved_quantity_key_, tag_inter_).UpdateDerivative(*S_, this->name(),
            key_, tag_inter_);
    db_->WriteVector("dTheta/du", S_->template GetDerivativePtr<CompositeVector>(conserved_quantity_key_, tag_inter_,
            key_, tag_inter_).ptr());

    f.Data()->reciprocal(S_->template GetDerivative<CompositeVector>(conserved_quantity_key_, tag_inter_, key_, tag_inter_));
    f.Data()->elementWiseMultiply(-1.0,
            S_->template Get<CompositeVector>(dudt_key_, tag_inter_), *f.Data(), 0.);
  }

 protected:
  using Base_t::plist_;
  using Base_t::domain_;
  using Base_t::tag_inter_;
  using Base_t::is_source_;
  using Base_t::S_;
  using Base_t::key_;
  using Base_t::mesh_;
  using Base_t::db_;
  using Base_t::vo_;
  using Base_t::conserved_quantity_key_;
  using Base_t::source_key_;

  Key dudt_key_;
  Key diffusion_key_;

};


// These are the actual classes.  They take the main PK
// (MixedFormParabolicPDE_*) and "mix in" the other classes that provide some
// functions.
//
// An implicitly time-integrated parabolic PDE in mixed form.
using PK_MixedFormParabolicPDE_Implicit =
    PK_Implicit_Adaptor<MixedFormParabolicPDE_Implicit<
                          PK_MixinConservationEquation<
                            PK_MixinImplicit<
                              PK_MixinLeafCompositeVector<
                                PK_Default>>>>>;

// An explicitly time-integrated parabolic PDE.  Maybe this should be renamed,
// it isn't really in mixed form anymore, but in primary form.
using PK_MixedFormParabolicPDE_Explicit =
    PK_Explicit_Adaptor<MixedFormParabolicPDE_Explicit<
                          PK_MixinConservationEquation<
                            PK_MixinExplicit<
                              PK_MixinLeafCompositeVector<
                                PK_Default>>>>>;

// A predictor-corrector scheme -- this combines both the implicit and explicit
// forms, then uses the explicit as a guess for the implicit.  The goal is for
// this to become the basis for higher order methods, but for now it really is
// just using the explicit as a guess for the implicit.
using PK_MixedFormParabolicPDE_PredictorCorrector =
    PK_ImplicitExplicit_Adaptor<MixedFormParabolicPDE_Implicit<
                          MixedFormParabolicPDE_Explicit<
                            PK_MixinConservationEquation<
                              PK_MixinPredictorCorrector<
                                PK_MixinLeafCompositeVector<
                                  PK_Default>>>>>>;


} // namespace Basic
} // namespace ATS
