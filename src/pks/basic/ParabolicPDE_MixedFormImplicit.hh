/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A parabolic PDE in mixed form for conservation, implicitly integrated.

/*!

This is the canonical nonlinear parabolic PDE, in mixed form.

.. math::
    \frac{\partial \Psi(u) }{\partial t} - \nabla \cdot K(u) \nabla \Phi(u) = Q(u,x,t)

where:

- :math:`u` the primary variable, key_
- :math:`\Psi` the conserved quantity, conserved_quantity_key_
- :math:`\Phi` the potential field, potential_key_
- :math:`K` the diffusion coefficient
- :math:`Q` any source term, source_key_

Any of these may be a function of the primary variable :math:`u`.

Note that frequently \Phi(u) = u (e.g. the potential field is the primary
variable) but not always -- specifically for surface water where the potential
field is :math:`\delta(p) + z` for primary variable pressure.


.._pk-parabolic-pde-mixed-form-implicit-spec
.. admonition:: pk-parabolic-pde-mixed-form-implicit-spec

   INCLUDES:
   - mixin-parabolic-pde-mixed-form-implicit-spec
   - mixin-conservation-equation-spec
   - mixin-implicit-spec
   - mixin-leaf-spec
   - pk-spec


.. _mixin-parabolic-pde-mixed-form-implicit-spec
.. admonition:: mixin-parabolic-pde-mixed-form-implicit-spec

    * `"potential key`" ``[string]`` The diffused potential field :math:`\Phi`

    * `"source key`" ``[string]`` **DOMAIN-source_sink** Units are in conserved
      quantity per second, :math:`Q`.

    * `"time discretization theta`" ``[double]`` **1.0** :math:`\theta` in a
      Crank-Nicholson time integration scheme.  1.0 implies fully implicit, 0.0
      implies explicit, 0.5 implies C-N.  Note, only used in the implicit
      scheme -- prefer to use an explicit time integrator over :math:`\theta ==
      0`.

    * `"modify predictor positivity preserving`" ``[bool]`` **false** If true,
      the predicted solution is modified to ensure that the conserved quantity
      is always > 0.  This does not guarantee positivity, so time integration
      method and/or timestep size may need to be chosen with care.

    * `"inverse`" ``[inverse-typed-spec]`` **optional**
      Linear inverse for the linear solve.  Only used if the time integration
      scheme is solved implicitly.


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
     <!-- this value is a small amount of water in mols/m^3 -->
     <Parameter name="absolute error tolerance" type="double" value="550.">
     <ParameterList name="inverse">
       ...
     </ParameterList>
     <ParameterList name="time integration">
       ...
     </ParameterList>
   </ParameterList>

*/

/*

DEVELOPER NOTES:

Note that PKs will mix this in with a time integration mixin and potentially
others to form a fully resolved PK.  For the simplest version of that, and for
examples of using this mixin, see PK_ParabolicPDE.hh

This Mixin fills parameters for a residual equation by combining three terms --
accumulation, diffusion, and (optionally) a source term.

- res_key_, the residual equation is an Evaluator_OperatorApply
- diffusion_key_, the div K grad term is an Evaluator_PDE_Diffusion of some form
- accumulation_key_, dPsi/dt above is an Evaluator_PDE_Accumulation, currently
  assumed to be lumped.
- source_key_, Q above

In the notation of Evaluator_OperatorApply, source_key_ and accumulation_key_
are both vectors, and are placed in the list of b vectors.

diffusion_key_ is an operator on the diagonal, e.g. A_00, and is applied to the
potential field Phi.

*/


#pragma once

#include "TreeVector.hh"
#include "Operator.hh"
#include "Operator_Factory.hh"
#include "Inverse.hh"
#include "SolverDefs.hh"

#include "PK_Adaptors.hh"
#include "PK_MixinImplicit.hh"
#include "PK_MixinLeaf.hh"
#include "PK_Default.hh"
#include "PK_Factory.hh"

#include "PK_MixinConservationEquation.hh"

namespace ATS {
namespace Basic {

using namespace Amanzi;

template <class Base_t>
class ParabolicPDE_MixedFormImplicit : public Base_t {

protected:
  using Base_t::plist_;
  using Base_t::domain_;
  using Base_t::tag_new_;
  using Base_t::tag_old_;
  using Base_t::S_;

  using Base_t::mesh_;
  using Base_t::db_;
  using Base_t::vo_;

  using Base_t::key_; // primary variable, u
  using Base_t::conserved_quantity_key_; // Psi
  using Base_t::is_source_;
  using Base_t::source_key_; // Q

  Key accumulation_key_; // dPsi/dt
  Key potential_key_; // Phi
  Key diffusion_key_; // div K grad operator
  Key res_key_; // residual

 public:
  // constructor from the base class is sufficient
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
    Teuchos::ParameterList& acc_list = S_->GetEvaluatorList(accumulation_key_);
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

    // Set up diffusion operator
    // -------------------------
    // This evaluator calculates q and local matrices for MFD or FV
    // representations of the div k K grad operator.  The evaluator itself is
    // the PDE_Diffusion(withGravity) object.
    //
    // Note we don't actually use the local matrices in this PK, so we don't
    // require them.  We require the global operator and that in turn requires
    // the local matrices.  But by setting up the parameter list here, we can
    // make the input file simpler and save the user from having to write all
    // this themselves.
    //
    // -- Input includes a potential variable, at the new time.
    potential_key_ = Keys::readKey(*plist_, domain_, "potential");
    if (potential_key_ != key_)
      S_->template RequireDerivative<CompositeVector,CompositeVectorSpace>(potential_key_,
              tag_new_, key_, tag_new_);
    // -- Require the Operator evaluator (a PDE_Diffusion object)
    diffusion_key_ = Keys::readKey(*plist_, domain_, "diffusion operator");
    Teuchos::ParameterList& diff_list = S_->GetEvaluatorList(diffusion_key_);
    if (!diff_list.isParameter("evaluator type")) diff_list.set("evaluator type", "diffusion operator");
    if (!diff_list.isParameter("local operator key")) diff_list.set("local operator key", diffusion_key_);
    if (!diff_list.isParameter("rhs key")) diff_list.set("rhs key", diffusion_key_+"_rhs");
    if (!diff_list.isParameter("flux key")) diff_list.set("flux key", diffusion_key_+"_flux");
    if (!diff_list.isParameter("boundary conditions key"))
      diff_list.set("boundary conditions key", diffusion_key_+"_bcs");
    diff_list.set("operator argument key", potential_key_);

    // Set up the residual evaluator
    // -----------------------------
    // Evaluator_OperatorApply sums the effects of operators, which are of the
    // form Ax = rhs, and vectors which are called b.
    //
    // In this case we have one operator, div K grad Phi, and two b vectors, -Q and dPsi/dt.
    //
    // See the documentation of Evaluator_OperatorApply for the details of these parameters.
    res_key_ = conserved_quantity_key_ + "_res";
    Teuchos::ParameterList& res_list = S_->GetEvaluatorList(res_key_);
    // -- residual is defined at the new time
    if (!res_list.isParameter("tag")) res_list.set("tag", tag_new_);
    // -- type is an Evaluator_OperatorApply, which uses an Operator object
    if (!res_list.isParameter("evaluator type")) res_list.set("evaluator type", "operator application");
    // -- the potential field is the diagonal primary variable, as this is a leaf PK
    res_list.set("diagonal primary x key", potential_key_);
    // -- the canonical nonlinear parabolic PDE has one operator, div K grad -
    if (!res_list.isParameter("diagonal local operators keys")) res_list.set("diagonal local operators keys",
                 Teuchos::Array<std::string>(1, diffusion_key_));
    // -- PDE_Diffusion stores its own RHS which must be added to the residual equation
    if (!res_list.isParameter("diagonal local operator rhs keys")) res_list.set("diagonal local operator rhs keys",
                 Teuchos::Array<std::string>(1, diffusion_key_+"_rhs"));

    // -- add in dPsi/dt and -Q
    Teuchos::Array<std::string> bs;
    Teuchos::Array<double> b_coefs;
    if (is_source_) {
      bs = std::vector<std::string>{ source_key_, accumulation_key_ };
      b_coefs = std::vector<double>{ -1.0, 1.0 };
    } else {
      bs = std::vector<std::string>{ accumulation_key_ };
      b_coefs = std::vector<double>{ 1.0 };
    }
    if (!res_list.isParameter("vector keys")) res_list.set("vector keys", bs);
    if (!res_list.isParameter("vector coefficients")) res_list.set("vector coefficients", b_coefs);

    // -- preconditioner -- could be provided in this PK's list, or in the
    //    State evaluator list.
    if (!res_list.isSublist("inverse")) res_list.set("inverse", plist_->sublist("inverse"));

    // NOTE: we cannot know the structure of the primary variable u, and
    // therefore the residual r here -- it may be CELL if FV, or it may be
    // CELL+FACE if MFD.  It will get set by the operator.  But we do need to
    // supply the mesh.
    S_->template Require<CompositeVector,CompositeVectorSpace>(key_, tag_new_)
        .SetMesh(mesh_);
    S_->template Require<CompositeVector,CompositeVectorSpace>(res_key_, tag_new_)
        .SetMesh(mesh_);
    S_->template RequireDerivative<Operators::Operator, Operators::Operator_Factory>(
      res_key_, tag_new_, potential_key_, tag_new_);

    // Now we can require the evaluator at the new time tag.
    S_->RequireEvaluator(res_key_, tag_new_);
  }


  void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                     Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> r)
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

    // evaluator r, assign it to the TreeVector data
    S_->GetEvaluator(res_key_, tag_new_).Update(*S_, this->name());
    r->Data()->assign(S_->template Get<CompositeVector>(res_key_, tag_new_));

    db_->WriteState(*S_, tag_old_);
    db_->WriteState(*S_, tag_new_);
    db_->WriteVector(res_key_, r->Data().ptr(), true);
  }

  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> r,
                          Teuchos::RCP<TreeVector> Pr)
  {
    Teuchos::OSTab tab = vo_->getOSTab();

    const auto& lin_op = S_->template GetDerivativePtr<Operators::Operator>(res_key_, tag_new_,
            potential_key_, tag_new_);

    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "Precon application:" << std::endl;
    db_->WriteVector("u_res", r->Data().ptr());
    int ierr = lin_op->applyInverse(*r->Data(), *Pr->Data());
    db_->WriteVector("PC*u_res", Pr->Data().ptr());

    if (potential_key_ != key_) {
      // u and Phi keys are NOT the same, must now apply the chain rule
      const auto& dpotential_dprimary = S_->template GetDerivative<CompositeVector>(potential_key_, tag_new_, key_, tag_new_);
      // NOTE: this could be lumped into a single kernel launch
      Pr->Data()->reciprocal(*Pr->Data());
      Pr->Data()->elementWiseMultiply(1.0, dpotential_dprimary, *Pr->Data(), 0.);
      Pr->Data()->reciprocal(*Pr->Data());
    }

    // return ierr;
    return 0;  // keep trying, even if linear solver fails to converge.
  }

  void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) {
    Teuchos::OSTab tab = vo_->getOSTab();
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "Precon update at t = " << t << std::endl;

    // evaluate dr/dPhi
    S_->GetEvaluator(res_key_, tag_new_).UpdateDerivative(*S_, this->name(), potential_key_, tag_new_);
    if (potential_key_ != key_)
      // evaluate dPhi/du
      S_->GetEvaluator(potential_key_, tag_new_).UpdateDerivative(*S_, this->name(), key_, tag_new_);
  }

  // -- Update diagnostics for vis.
  void CalculateDiagnostics(const Key& tag) {}

};


/*

DEVELOPER NOTE:

This is the actual PK class.  They take the main PK
(ParabolicPDE_MixedFormImplicit) and "mix in" the other classes that provide
capability for:

- PK_MixinConservationEquation: a conserved quantity and an error norm
- PK_MixinImplicit: implicit time integration
- PK_MixinLeaf: a domain and mesh, debugger
- PK_Default: a verbose object and name

*/
using PK_ParabolicPDE_MixedFormImplicit =
    PK_Implicit_Adaptor<ParabolicPDE_MixedFormImplicit<
                          PK_MixinConservationEquation<
                            PK_MixinImplicit<
                              PK_MixinLeafCompositeVector<
                                PK_Default>>>>>;

} // namespace Basic
} // namespace ATS
