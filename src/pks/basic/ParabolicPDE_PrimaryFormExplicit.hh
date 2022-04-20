/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A parabolic PDE in primary form, explicitly integrated.

/*!

This is the canonical nonlinear parabolic PDE, in primary form for explicit
time integration.

.. math::
    \frac{d \Psi(u) }{d u} \frac{\partial u}{\partial t} - \nabla \cdot K(u) \nabla u = Q(u,x,t)

where:

- :math:`u` the primary variable, key_
- :math:`\Phi` the potential field, potential_key_
- :math:`K` the diffusion coefficient
- :math:`Q` any source term, source_key_
o
Note that this, unlike ParabolicPDE_MixedFormImplicit, does not currently
support a potential field that is not u.

Any of these may be a function of the primary variable :math:`u`.


.._pk-parabolic-pde-primary-form-explicit-spec
.. admonition:: pk-parabolic-pde-primary-form-explicit-spec

   INCLUDES:
   - mixin-parabolic-pde-primary-form-explicit-spec
   - mixin-conservation-equation-spec
   - mixin-explicit-spec
   - mixin-leaf-spec
   - pk-spec


.. _mixin-parabolic-pde-primary-form-explicit-spec
.. admonition:: mixin-parabolic-pde-primary-form-explicit-spec

    * `"source key`" ``[string]`` **DOMAIN-source_sink** Units are in conserved
      quantity per second, :math:`Q`.


Example: Richards equation
~~~~~~~~~~~~~~~~~~~~~~~~~~

In the case of Richards equation, the user would likely set the following
names as parameters:

.. xml::
   <ParameterList name="Richards">
     <Parameter name="domain" type="string" value="domain">
     <Parameter name="primary variable key" type="string" value="pressure">
     <Parameter name="conserved quantity key" type="string" value="water_content">
     <Parameter name="potential key" type="string" value="pressure">
     <Parameter name="source key" type="string" value="water_source_sink">
     <ParameterList name="time integration">
       ...
     </ParameterList>
   </ParameterList>

*/

#pragma once

#include "TreeVector.hh"
#include "Operator.hh"
#include "Operator_Factory.hh"
#include "Inverse.hh"
#include "SolverDefs.hh"

#include "PK_Adaptors.hh"
#include "PK_MixinExplicit.hh"
#include "PK_MixinLeaf.hh"
#include "PK_Default.hh"
#include "PK_Factory.hh"

#include "PK_MixinConservationEquation.hh"

namespace ATS {
namespace Basic {

using namespace Amanzi;

template <class Base_t>
class ParabolicPDE_PrimaryFormExplicit : public Base_t {

 protected:
  using Base_t::plist_;
  using Base_t::domain_;
  using Base_t::tag_inter_; // tag at which the RHS of dudt = ... is evaluated
  using Base_t::is_source_;
  using Base_t::S_;

  using Base_t::mesh_;
  using Base_t::db_;
  using Base_t::vo_;

  using Base_t::key_; // primary variable, u
  using Base_t::conserved_quantity_key_; // Psi
  using Base_t::source_key_; // Q

  Key dudt_key_; // time derivative of the primary variable
  Key diffusion_key_; // div K grad operator

 public:
  // base constructor is sufficient
  using Base_t::Base_t;

  void Setup()
  {
    Base_t::Setup();
    Base_t::SetupAtTag(tag_inter_);

    // set up the diffusion operator
    diffusion_key_ = Keys::readKey(*plist_, domain_, "diffusion operator");
    Teuchos::ParameterList& diff_list = S_->GetEvaluatorList(diffusion_key_);
    if (!diff_list.isParameter("evaluator type")) diff_list.set("evaluator type", "diffusion operator");
    if (!diff_list.isParameter("local operator key")) diff_list.set("local operator key", diffusion_key_);
    if (!diff_list.isParameter("rhs key")) diff_list.set("rhs key", diffusion_key_+"_rhs");
    if (!diff_list.isParameter("boundary conditions key"))
      diff_list.set("boundary conditions key", diffusion_key_+"_bcs");
    diff_list.set("operator argument key", key_);

    // require Psi at the intermediate/old time
    S_->template Require<CompositeVector,CompositeVectorSpace>(conserved_quantity_key_, tag_inter_)
        .SetMesh(mesh_)->AddComponent("cell", Amanzi::AmanziMesh::CELL, 1);
    // require dPsi/du at the intermediate/old time
    S_->template RequireDerivative<CompositeVector,CompositeVectorSpace>(conserved_quantity_key_, tag_inter_,
            key_, tag_inter_);
    S_->RequireEvaluator(conserved_quantity_key_, tag_inter_);

    // require du/dt
    dudt_key_ = key_ + "_t";
    // ETC: this used to be here, likely a bug?
    //dudt_key_ = conserved_quantity_key_ + "_t";
    Teuchos::ParameterList& res_list = S_->GetEvaluatorList(dudt_key_);
    // -- operator
    if (!res_list.isParameter("tag")) res_list.set("tag", tag_inter_);
    if (!res_list.isParameter("evaluator type")) res_list.set("evaluator type", "operator application");
    if (!res_list.isParameter("diagonal primary x key")) res_list.set("diagonal primary x key", key_);
    if (!res_list.isParameter("diagonal local operators keys")) res_list.set("diagonal local operators keys",
                 Teuchos::Array<std::string>(1, diffusion_key_));
    if (!res_list.isParameter("diagonal local operator rhs keys")) res_list.set("diagonal local operator rhs keys",
                 Teuchos::Array<std::string>(1, diffusion_key_+"_rhs"));

    // -- rhs for source and accumulation term
    Teuchos::Array<std::string> bs;
    Teuchos::Array<double> b_coefs;
    if (is_source_) {
      bs = std::vector<std::string>{ source_key_ };
      b_coefs = std::vector<double>{ -1.0 };
    }
    if (!res_list.isParameter("vector keys")) res_list.set("vector keys", bs);
    if (!res_list.isParameter("vector coefficients")) res_list.set("vector coefficients", b_coefs);

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

};

/*

DEVELOPER NOTE:

This is the actual PK class.  They take the main PK
(ParabolicPDE_MixedFormImplicit) and "mix in" the other classes that provide
capability for:

- PK_MixinConservationEquation: a conserved quantity and an error norm
- PK_MixinExplicit: explicit time integration
- PK_MixinLeaf: a domain and mesh, debugger
- PK_Default: a verbose object and name

*/
using PK_ParabolicPDE_PrimaryFormExplicit =
    PK_Explicit_Adaptor<ParabolicPDE_PrimaryFormExplicit<
                          PK_MixinConservationEquation<
                            PK_MixinExplicit<
                              PK_MixinLeafCompositeVector<
                                PK_Default>>>>>;

} // namespace Basic
} // namespace ATS
