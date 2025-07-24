/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
/*!

This MPC implements an operator-split coupled water, which splits the lateral
overland fluxes from surface sources and subsurface flow.

To advance the coupled flow problem from :math:`t^n` to :math:`t^{n + 1}`, we
perform two steps.  First, solve the diffusion wave equation:

.. math::

   \frac{\partial \Theta_s(p_s^*)}{\partial t} = \nabla k_s \cdot \nabla (z+h(p_s^*))

for surface pressure :math:`p_s^*`.  This is colloquially called the "star"
system.  Then, one of three algorithms is performed for the "primary" system,
which is a variation on :ref:`Integrated Hydrology` without lateral surface
fluxes and with extra coupling to the star system.

**Pressure Coupling**

`"pressure`" coupling is a standard operator splitting.  :math:`p_s^*` becomes
the initial value for the coupled surface sources and subsurface equations:

.. math::

   \frac{\partial \Theta_s(p_s)}{\partial t} &= Q_ext + q_{ss} \\
   \frac{\partial \Theta(p)}{\partial t} &= \nabla \cdot k K (\nabla p + \rho g \hat{z}) \\
   -k K (\nabla p + \rho g \hat{z}) \cdot \hat{n} |_\Gamma &= q_{ss} \\
   p|_\Gamma &= p_s \\
   p_s(t^n) &= p_s^*

Note the lack of a lateral flow term in the overland equation (relative to the
standard diffusion wave equation shown in `Integrated Hydrology`_).  This
system is coupled and solved in the same discrete way as `Integrated
Hydrology`_.

**Flux Coupling**

`"flux`" coupling, rather than setting the initial pressure from the
:math:`p_s^*`, instead provides the divergence of fluxes in the lateral flow
system as a source to the surface system:

.. math::

   \frac{\partial \Theta_s(p_s)}{\partial t} &= \frac{\partial \Theta_s(p_s^*)}{\partial t} + Q_ext + q_{ss} \\
   \frac{\partial \Theta(p)}{\partial t} &= \nabla \cdot k K (\nabla p + \rho g \hat{z}) \\
   -k K (\nabla p + \rho g \hat{z}) \cdot \hat{n} |_\Gamma &= q_{ss} \\
   p|_\Gamma &= p_s

The advantage of this approach is that it more stably handles the case of
wetting up -- when a dry surface cell first gets overland flow into that cell,
it requires that :math:`p_s > p_{atm}`.  But if the subsurface below it is
unsaturated, this can create a large gradient in pressure that will immediately
be eliminated in the subsurface solve once "run-on" infiltrates.  This jump
between a :math:`p_s^* > p_atm` and :math:`p_s < p_atm` is unstable, and hard
on the nonlinear solver.

Imposing the run-on as a source of water rather than as a initial pressure is
much more stable for run-on.


**Hybrid Coupling**

While the `"flux`" approach is more stable than `'pressure`" in cells
experiencing run-on, it is not necessary, and potentially problematic in the
case of run-off.  In those cases the `"pressure`" case is more stable.
Therefore, a `"hybrid`" coupling approach is used most frequently; this uses
the `"pressure`' algorithm where the divergence of lateral surface fluxes is
negative (e.g. run-off) and the `"flux`" algorithm elsewhere.

All three algorithms should result in the same (converged) solution.  The
`"hybrid`" algorithm is the most robust for numerical performance.


Additionally, the subsurface domain may be treated as either a 3D domain
(solving 3D Richards equations) or as a domain-set of many 1D, vertical
columns.  In this case, the second system of equations is implemented on each
subsurface column individually -- there is nothing coupling these columns in
the second system.  Lateral subsurface flow is ignored in this case.  This
allows subcycling of individual columns, and a much more efficient solve.  This
is most appropriate at larger or "intermediate" scales, where lateral flow in
the subsurface is small.

`"PK type`" = `"operator split coupled water`"

.. _pk-operator-split-coupled-water-spec:
.. admonition:: pk-operator-split-coupled-water-spec

   * `"PKs order`" ``[Array(string)]`` Order is {star_system, primary_system}.
     Note that the sub-PKs are likely a :ref:`Overland Flow PK` for the "star"
     system and a :ref:`Integrated Hydrology` MPC for the "primary" system.

   * `"domain name`" ``[string]`` The subsurface domain, e.g. `"domain`" (for a
     3D subsurface) or `"column:*`" (for the intermediate scale, columnar model).

   * `"star domain name`" ``[string]`` The surface domain, typically
     `"surface_star`" by convention.

   * `"coupling type`" ``[string]`` **"hybrid"** One of: `"pressure`", `"flux`",
     or `"hybrid`" (see above).

   INCLUDES:

   - ``[mpc-spec]`` *Is an* :ref:`MPC`.
   - ``[mpc-subcycled-spec]`` *Is a* :ref:`Subcycling MPC`.

*/

#pragma once

#include "PK.hh"
#include "mpc_subcycled.hh"

namespace Amanzi {

class MPCCoupledWaterSplitFlux : public MPCSubcycled {
 public:
  MPCCoupledWaterSplitFlux(Teuchos::ParameterList& FElist,
                           const Teuchos::RCP<Teuchos::ParameterList>& plist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution);

  // PK methods
  void parseParameterList() override;

  // -- initialize in reverse order
  virtual void Initialize() override;
  virtual void Setup() override;

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;

 protected:
  void CopyPrimaryToStar_();
  void CopyStarToPrimary_();

  void CopyPrimaryToStar_DomainSet_();
  void CopyStarToPrimary_DomainSet_Pressure_();
  void CopyStarToPrimary_DomainSet_Flux_();
  void CopyStarToPrimary_DomainSet_Hybrid_();

  void CopyPrimaryToStar_Standard_();
  void CopyStarToPrimary_Standard_Pressure_();
  void CopyStarToPrimary_Standard_Flux_();
  void CopyStarToPrimary_Standard_Hybrid_();

  Tag get_ds_tag_next_(const std::string& subdomain)
  {
    if (subcycling_[1])
      return Tag{ Keys::cleanName(tags_[1].second.get() + "_" +
                                  Keys::getDomainSetIndex(subdomain)) };
    else return Tag{ tags_[1].second };
  }
  Tag get_ds_tag_current_(const std::string& subdomain)
  {
    if (subcycling_[1])
      return Tag{ Keys::cleanName(tags_[1].first.get() + "_" +
                                  Keys::getDomainSetIndex(subdomain)) };
    else return Tag{ tags_[1].first };
  }


 protected:
  Key p_primary_variable_;
  Key p_primary_variable_suffix_;
  Key p_sub_primary_variable_;
  Key p_sub_primary_variable_suffix_;
  Key p_primary_variable_star_;
  Key p_conserved_variable_;
  Key p_conserved_variable_star_;
  Key p_lateral_flow_source_;
  Key p_lateral_flow_source_suffix_;

  Key cv_key_;

  Key domain_set_;
  Key domain_;
  Key domain_sub_;
  Key domain_star_;

  std::string coupling_;

  bool is_domain_set_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCCoupledWaterSplitFlux> reg_;
};

} // namespace Amanzi
