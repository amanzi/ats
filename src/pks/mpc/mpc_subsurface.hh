/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A coupler which solves flow and energy in the subsurface.
/*!

This MPC provides most nearly all terms for an approximate Jacobian for
coupling three-phase Richards equation (the `Permafrost Flow PK`_) to the
three-phase Energy equation (the `Three-Phase subsurface Energy PK`_).

Many options are provided for turning on and off various aspects of this
Jacobian, so it is useful to mathematically write out these terms.  The
equations are:

.. math::
    \frac{\partial \Theta}{\partial t} - \nabla \frac{k_r n_l}{\mu} K ( \nabla p + \rho g \cdot \hat{z} ) = Q_w \\
    \frac{\partial E}{\partial t} - \nabla \cdot \kappa \nabla T + \nabla \cdot \mathbf{q} e(T) = Q_w e(T) + Q_e

Note that all of the following are dependent on :math:`p` and/or :math:`T`:

.. math::
    \Theta(p,T), k_r(p,T), n_l(p,T), \mu(T), \rho(p,T), E(p,T), \kappa(p,T), e(T)

Also, both source terms :math:`Q_w` and :math:`Q_e` may or may not depend on :math:`p` and :math:`T`.

Note also that the Darcy flux :math:`\mathbf{q}` used in the advection of energy is given by the Darcy flux:

.. math::
    \mathbf{q} = -\frac{k_r n_l}{\mu} K ( \nabla p + \rho g \cdot \rho g \cdot \hat{z} )

Differentiating these two equations in their two unknowns gives the following four blocks in the approximate Jacobian:

:math:`\frac{\partial F_1}{\partial p}`: this is the Richards equation diagonal block, and is controlled inside that PK.

:math:`\frac{\partial F_1}{\partial T}` includes terms for:

- :math:`\frac{\partial \Theta}{\partial T}` This term is the cell-local diagonal block.

- The partial derivative of the divergence of the Darcy flux with respect to
  temperature is dominated by :math:`\frac{\partial}{\partial T} \frac{k_r
  n_l}{\mu}`.  This is because the relative permeability is strongly
  dependent upon phase change (the freezing equals drying approximation).  This
  term is referred to as the "d div q / dT" term.

:math:`\frac{\partial F_2}{\partial p}` includes terms for:

- :math:`\frac{\partial E}{\partial p}` This term is the cell-local diagonal block.

- The partial derivative of the energy diffusion term with respect to pressure
  is dominated by :math:`\frac{\partial \kappa}{\partial p}` through phase
  change -- at a constant temperature, but changing pressure, phase change can
  result in large changes to thermal conductivity.  This is referred to as the
  "div K grad T / dp" term.

:math:`\frac{\partial F_2}{\partial T}`: this is the energy equation diagonal block, and is controlled inside that PK.

Also, at this level, where we know more about the flux used in the energy
equation (it is the Darcy flux), we can do a better approximation of the
derivative of the advection of energy term with respect to both temperature and
pressure.  For instance, enthalpy is only weakly dependent on pressure, so we
can use the derivative of the divergence of the Darcy flux with respect to
pressure (from the Richards block) in the advection term in the
:math:`\frac{\partial F_2}{\partial p}` block, and approximate
:math:`\frac{\partial k_r}{\partial T}` in the advection term as well.  These
terms are referred to as "div hq / dp,T terms".  Note the missing initial "d"
here relative to other terms.

The behavior of this MPC's preconditioner can be set by an option,
`"preconditioner type`".  Really users should not change this from the default,
except in expert cases or for comparison's sake, but the options are:

- `"picard`" is the default, this uses all available terms, and enables the
  "suppress" options for finer-grained control.

- `"none`" No preconditioner never works.

- `"block diagonal`" This is what one would get from the default StrongMPC_.  This probably never works.

- `"no flow coupling`" This keeps the accumulation terms, but turns off all the
  non-local blocks.  This is equivalent to `Coupled Cells MPC`_.

- `"ewc`" **CURRENTLY DEPRECATED/BROKEN/DISABLED** In addition to the
  `"picard`" coupling, this also *always* does a change of variables, whereby
  we first invert to calculate primary variable corrections, then do a change
  of variables to calculate the linearized corrections in energy and water
  content space.  We then apply those corrections, and invert to find the
  primary variable changes that would have made those corrections.  This is
  called the "energy and water content" algorithm, and is related to similar
  variable changing approaches by Krabbenhoft (for flow) and Knoll (for
  energy), but in the multivariate approach.  This is somewhat bad, becuase
  while it fixes some corrections, it breaks others.

- `"smart ewc`" **CURRENTLY DEPRECATED/BROKEN/DISABLED** Does the `"ewc`"
  algorithm above, but tries to be smart about when to do it.  This algorithm
  helps when we are about to fall off of the latent heat cliff.  If we can
  guess when to do it, we have a better chance of not breaking things.  This
  seems like it ought to be helpful, but often doesn't do as much as one might
  hope.


Note this "ewc" algorithm is just as valid, and more useful, in the predictor
(where it is not deprecated/disabled).  There, we extrapolate a change in
pressure and temperature, but often do better to extrapolate in water content
and energy space, then invert (locally) for pressure and temperature
corrections that meet that extrapolation.  Both of these globalization
algorithms are supported by the `EWC Globalization Delegate`_ object.

.. _mpc-subsurface-spec:
.. admonition:: mpc-subsurface-spec

    * `"domain name`" ``[string]`` Domain of simulation

    * `"preconditioner type`" ``[string]`` **picard** See the above for
      detailed descriptions of the choices.  One of: `"none`", `"block
      diagonal`", `"no flow coupling`", `"picard`", `"ewc`", and `"smart ewc`".

    * `"supress Jacobian terms: div hq / dp,T`" ``[bool]`` **false** If using picard or ewc, do not include this block in the preconditioner.
    * `"supress Jacobian terms: d div q / dT`" ``[bool]`` **false** If using picard or ewc, do not include this block in the preconditioner.
    * `"supress Jacobian terms: d div K grad T / dp`" ``[bool]`` **false** If using picard or ewc, do not include this block in the preconditioner.

    * `"ewc delegate`" ``[mpc-delegate-ewc-spec]`` A `EWC Globalization Delegate`_ spec.

    INCLUDES:

    - ``[strong-mpc-spec]`` *Is a* StrongMPC_.

 */

#ifndef MPC_SUBSURFACE_HH_
#define MPC_SUBSURFACE_HH_

#include "TreeOperator.hh"
#include "pk_physical_bdf_default.hh"
#include "strong_mpc.hh"

namespace Amanzi {

class MPCDelegateEWCSubsurface;

namespace Operators {
class PDE_Diffusion;
class PDE_DiffusionWithGravity;
class PDE_Accumulation;
class Operator;
class UpwindTotalFlux;
class UpwindArithmeticMean;
class Upwinding;
} // namespace Operators

namespace Flow {
class Richards;
}

class MPCSubsurface : public StrongMPC<PK_PhysicalBDF_Default> {
 public:
  MPCSubsurface(Teuchos::ParameterList& pk_tree_list,
                const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& soln);

  // -- Initialize owned (dependent) variables.
  virtual void Setup() override;
  virtual void Initialize() override;
  virtual void set_tags(const Tag& tag_current, const Tag& tag_next) override;

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h,
                               Teuchos::RCP<const TreeVector> up0,
                               Teuchos::RCP<TreeVector> up) override;

  // -- Update the preconditioner to be physically consistent
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

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
  Teuchos::RCP<Operators::PDE_DiffusionWithGravity> ddivq_dT_;
  Teuchos::RCP<Operators::Upwinding> upwinding_dkrdT_;
  // -- d ( dWC/dt ) / dT terms
  Teuchos::RCP<Operators::PDE_Accumulation> dWC_dT_;

  // dE / dp off-diagonal block
  Teuchos::RCP<Operators::Operator> dE_dp_block_;
  // -- d ( div K grad T ) / dp terms
  Teuchos::RCP<Operators::PDE_Diffusion> ddivKgT_dp_;
  Teuchos::RCP<Operators::Upwinding> upwinding_dkappa_dp_;
  // -- d ( div hq ) / dp terms
  Teuchos::RCP<Operators::PDE_DiffusionWithGravity> ddivhq_dp_;
  Teuchos::RCP<Operators::UpwindTotalFlux> upwinding_hkr_;
  Teuchos::RCP<Operators::UpwindTotalFlux> upwinding_dhkr_dp_;
  // -- d ( dE/dt ) / dp terms
  Teuchos::RCP<Operators::PDE_Accumulation> dE_dp_;

  // dE / dT on-diagonal block additional terms that use q info
  // -- d ( div hq ) / dT terms
  Teuchos::RCP<Operators::PDE_DiffusionWithGravity> ddivhq_dT_;
  Teuchos::RCP<Operators::UpwindTotalFlux> upwinding_dhkr_dT_;

  // friend sub-pk Richards (need K_, some flags from private data)
  //Teuchos::RCP<Flow::Richards> richards_pk_;

  Key domain_name_;
  Key temp_key_;
  Key pres_key_;
  Key e_key_;
  Key wc_key_;
  Key tc_key_;
  Key uw_tc_key_;
  Key kr_key_;
  Key uw_kr_key_;
  Key enth_key_;
  Key hkr_key_;
  Key uw_hkr_key_;
  Key energy_flux_key_;
  Key water_flux_key_;
  Key water_flux_dir_key_;
  Key rho_key_;
  Key duw_krdT_key_;
  Key duw_tcdp_key_;

  bool is_fv_;

  // EWC delegate
  Teuchos::RCP<MPCDelegateEWCSubsurface> ewc_;

  // cruft for easier global debugging
  bool dump_;
  int update_pcs_;
  Teuchos::RCP<Debugger> db_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCSubsurface> reg_;
};

} // namespace Amanzi

#endif
