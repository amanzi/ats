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
    \frac{\partial \Phi(u) }{\partial t} = Q(u,x,t)

Where the conserved quantity :math:`\Phi` is a function of the primary variable
:math:`u`, and a source term :math:`Q` is provided.
    
.. _conservation-ode-pk-spec:
.. admonition:: conservation-ode-pk-spec

    * `"domain`" ``[string]`` Mesh on which the balance is to be done.

    * `"primary variable key`" ``[string]`` The primary variable, :math:`u`.
      Note there is no default -- this must be provided by the user.

    * `"conserved quantity key`" ``[string]`` The conserved quantity :math:`\Phi`

    * `"source key`" ``[string]`` **DOMAIN-source_sink** Units are in conserved
      quantity per second, :math:`Q`.

    * `"time discretization theta`" ``[double]`` **1.0** :math:`\theta` in a
      Crank-Nicholson time integration scheme.  1.0 implies fully implicit, 0.0
      implies explicit, 0.5 implies C-N.  Note, only used in the implicit
      scheme.

    * `"modify predictor positivity preserving`" ``[bool]`` **false** If true,
      predictors are modified to ensure that the conserved quantity is always >
      0.  These systems may be stiff, and this does not guarantee positivity,
      so time integration methods may need to be chosen with care.

    * `"absolute error tolerance`" ``[double]`` **550.0** a_tol in the standard
      error norm calculation.  Defaults to a small amount of water.  Units are
      the same as the conserved quantity.

    * `"preconditioner`" ``[preconditioner-typed-spec]`` **optional**
      Preconditioner for the solve.  Only used if the time integration scheme
      is solved implicitly.

    * `"linear solver`" ``[linear-solver-typed-spec]`` **optional** May be used
      to improve the inverse of the preconditioner.  Only used if this PK is
      implicitly solved, but not implicitly coupled at a higher level.  See
      LinearOperator_.

*/

#pragma once

#include "TreeVector.hh"
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
class ConservationODE_Implicit : public Base_t {

public:
  using Base_t::Base_t;

  void Setup() {
    Base_t::Setup();
    Base_t::SetupAtTag(tag_new_);

    // set up the accumulation evaluator
    accumulation_key_ = conserved_quantity_key_ + "_t";
    Teuchos::ParameterList& acc_list = S_->FEList().sublist(accumulation_key_);
    acc_list.set("conserved quantity key", conserved_quantity_key_);
    acc_list.set("tag", tag_new_);
    acc_list.set("tag old", tag_old_);
    acc_list.set("tag new", tag_new_);
    acc_list.set("evaluator type", "accumulation operator");

    S_->template Require<CompositeVector,CompositeVectorSpace>(accumulation_key_, tag_new_)
        .SetMesh(mesh_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
    S_->RequireEvaluator(accumulation_key_, tag_new_);

    // set up the residual evaluator
    res_key_ = key_ + "_res";
    Teuchos::ParameterList& res_list = S_->FEList().sublist(res_key_);
    res_list.set("tag", tag_new_);
    res_list.set("evaluator type", "additive");
    std::vector<std::string> deps{accumulation_key_, source_key_};
    res_list.set("dependencies", Teuchos::Array<std::string>(deps));
    res_list.set(accumulation_key_+" coefficient", 1.0);
    res_list.set(source_key_+" coefficient", -1.0);

    S_->template Require<CompositeVector,CompositeVectorSpace>(res_key_, tag_new_)
        .SetMesh(mesh_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
    S_->template RequireDerivative<CompositeVector,CompositeVectorSpace>(res_key_, tag_new_,
            key_, tag_new_);
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
    db_->WriteVectors({"u_old","u_new"}, {u_old->Data().ptr(), u_new->Data().ptr()});
    
    S_->GetEvaluator(res_key_, tag_new_).Update(*S_, this->name());
    f->Data()->assign(S_->template Get<CompositeVector>(res_key_, tag_new_));
    db_->WriteVector("u_res", f->Data().ptr());
  }

  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> r,
                          Teuchos::RCP<TreeVector> Pr)
  {
    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "Precon application:" << std::endl;
    
    db_->WriteVector("u_res", r->Data().ptr());

    const auto& dr_du = S_->template GetDerivative<CompositeVector>(res_key_, tag_new_,
            key_, tag_new_);
    Pr->Data()->reciprocal(dr_du);
    Pr->Data()->elementWiseMultiply(1.0, *r->Data(), *Pr->Data(), 0.);

    db_->WriteVector("PC*u_res", Pr->Data().ptr());
    return 0;
  }

  void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "Precon update at t = " << t << std::endl;
    
    S_->GetEvaluator(res_key_, tag_new_).UpdateDerivative(*S_, this->name(),
            key_, tag_new_);

    db_->WriteVector("  diagonal PC",
                     S_->template GetDerivative<CompositeVector>(res_key_, tag_new_,
                             key_, tag_new_));
  }

  // -- Update diagnostics for vis.
  void CalculateDiagnostics(const Key& tag) {}

 protected:
  using Base_t::tag_new_;
  using Base_t::tag_old_;
  using Base_t::S_;
  using Base_t::key_;
  using Base_t::mesh_;
  using Base_t::db_;
  using Base_t::vo_;
  using Base_t::conserved_quantity_key_;
  using Base_t::source_key_;

  Key accumulation_key_;
  Key res_key_;

};


template <class Base_t>
class ConservationODE_Explicit : public Base_t {

public:
  using Base_t::Base_t;

  void Setup() {
    Base_t::Setup();
    Base_t::SetupAtTag(tag_inter_);
  
    // here we just need the source and the conserved quantity's derivative
    S_->template Require<CompositeVector,CompositeVectorSpace>(source_key_, tag_inter_)
        .SetMesh(mesh_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
    S_->RequireEvaluator(source_key_, tag_inter_);

    S_->template Require<CompositeVector,CompositeVectorSpace>(conserved_quantity_key_, tag_inter_)
        .SetMesh(mesh_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
    S_->template RequireDerivative<CompositeVector,CompositeVectorSpace>(conserved_quantity_key_, tag_inter_,
            key_, tag_inter_);
    S_->RequireEvaluator(conserved_quantity_key_, tag_inter_);
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

    S_->GetEvaluator(source_key_, tag_inter_).Update(*S_, this->name());
    db_->WriteVector("source", S_->template GetPtr<CompositeVector>(source_key_, tag_inter_).ptr());

    S_->GetEvaluator(conserved_quantity_key_, tag_inter_).UpdateDerivative(*S_, this->name(),
            key_, tag_inter_);
    db_->WriteVector("dTheta/du", S_->template GetDerivativePtr<CompositeVector>(conserved_quantity_key_, tag_inter_,
            key_, tag_inter_).ptr());

    f.Data()->reciprocal(S_->template GetDerivative<CompositeVector>(conserved_quantity_key_, tag_inter_, key_, tag_inter_));
    f.Data()->elementWiseMultiply(1.0, S_->template Get<CompositeVector>(source_key_, tag_inter_), *f.Data(), 0.);
  }


 protected:
  using Base_t::tag_inter_;
  using Base_t::S_;
  using Base_t::key_;
  using Base_t::mesh_;
  using Base_t::db_;
  using Base_t::vo_;
  using Base_t::conserved_quantity_key_;
  using Base_t::source_key_;

};

  




using PK_ConservationODE_Implicit =
    PK_Implicit_Adaptor<ConservationODE_Implicit<
                          PK_MixinConservationEquation<
                            PK_MixinImplicit<
                              PK_MixinLeafCompositeVector<
                                PK_Default>>>>>;

using PK_ConservationODE_Explicit =
    PK_Explicit_Adaptor<ConservationODE_Explicit<
                          PK_MixinConservationEquation<
                            PK_MixinExplicit<
                              PK_MixinLeafCompositeVector<
                                PK_Default>>>>>;

using PK_ConservationODE_PredictorCorrector =
    PK_ImplicitExplicit_Adaptor<ConservationODE_Implicit<
                          ConservationODE_Explicit<
                            PK_MixinConservationEquation<
                              PK_MixinPredictorCorrector<
                                PK_MixinLeafCompositeVector<
                                  PK_Default>>>>>>;


} // namespace Basic
} // namespace ATS

