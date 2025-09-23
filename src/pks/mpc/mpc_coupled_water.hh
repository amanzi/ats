/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

     \frac{\partial \Theta_s(p_s^*)}{\partial t} &= \nabla k_s \cdot \nabla (z+h(p_s^*)) + Q_{ext} + q_{ss} \\
     \frac{\partial \Theta(p)}{\partial t} &= \nabla \cdot k K (\nabla p + \rho g \hat{z}) \\
     & -k K (\nabla p + \rho g \hat{z}) \cdot \hat{n} |_{\Gamma} = q_{ss} \\
     & p|_{\Gamma} = p_s
*/
/*!

This MPC couples Richards equation (subsurface flow) to the diffusion wave
equation (surface flow) through continuity of both pressure and fluxes.  In its
standard form, it couples pressure at a the surface of the domain,
:math:`\Gamma`, by solving:

.. math::

   \frac{\partial \Theta_s(p_s^*)}{\partial t} &= \nabla k_s \cdot \nabla (z+h(p_s^*)) + Q_{ext} + q_{ss} \\
   \frac{\partial \Theta(p)}{\partial t} &= \nabla \cdot k K (\nabla p + \rho g \hat{z}) \\
   -k K (\nabla p + \rho g \hat{z}) \cdot \hat{n} |_{\Gamma} &= q_{ss} \\
   p|_{\Gamma} &= p_s

:math:`q_{ss}`, the "exfiltration" of water from subsurface to surface, is
eliminated discretely (enforcing flux continuity) and the pressure equality
condition (the last equation) is enforced by elimintating the surface unknowns
discretely.

This leverages subsurface discretizations that include face-based unknowns, and
notes that those face unknowns that correspond to surface faces are co-located
with the surface cell pressure, and therefore are equivalent.  In this approach
(see `Coon et al WRR 2020 <https://doi.org/10.1016/j.advwatres.2020.103701>`_),
the surface equations are directly assembled into the subsurface discrete
operator.

As a result, this requires a subsurface discretization that uses face-based
unknowns, e.g. any of the `"mfd: *`" class of methods.

`"PK type`" = `"coupled water`"

.. _pk-coupled-water-spec:
.. admonition:: pk-coupled-water-spec

   * `"PKs order`" ``[Array(string)]`` Order must be {subsurface_flow_pk,
     surface_flow_pk}.

   * `"subsurface domain name`" ``[string]`` **domain**

   * `"surface domain name`" ``[string]`` **surface**

   * `"water delegate`" ``[mpc-delegate-water-spec]`` A :ref:`Coupled Water
     Globalization Delegate` spec.

   INCLUDES:

   - ``[strong-mpc-spec]`` *Is a* :ref:`Strong MPC`.

*/

#ifndef PKS_MPC_COUPLED_WATER_HH_
#define PKS_MPC_COUPLED_WATER_HH_

#include "Operator.hh"
#include "mpc_delegate_water.hh"
#include "pk_physical_bdf_default.hh"

#include "strong_mpc.hh"

namespace Amanzi {

class MPCCoupledWater : public StrongMPC<ATS_Physics::PK_PhysicalBDF_Default> {
 public:
  MPCCoupledWater(Teuchos::ParameterList& FElist,
                  const Teuchos::RCP<Teuchos::ParameterList>& plist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  // Parse the local parameter list and add entries to the global list
  virtual void parseParameterList() override;

  virtual void Setup() override;
  virtual void Initialize() override;

  // -- computes the non-linear functional g = g(t,u,udot)
  //    By default this just calls each sub pk FunctionalResidual().
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g) override;

  // -- Apply preconditioner to u and returns the result in Pu.
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> Pu) override;

  // -- Modify the predictor.
  virtual bool ModifyPredictor(double h,
                               Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) override;

  // -- Modify the correction.
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult ModifyCorrection(
    double h,
    Teuchos::RCP<const TreeVector> res,
    Teuchos::RCP<const TreeVector> u,
    Teuchos::RCP<TreeVector> du) override;

  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> res) override;

  Teuchos::RCP<Operators::Operator> preconditioner() { return precon_; }

 protected:
  // void
  // UpdateConsistentFaceCorrectionWater_(const Teuchos::RCP<const TreeVector>& u,
  //         const Teuchos::RCP<TreeVector>& Pu);

 protected:
  std::string domain_surf_, domain_ss_;
  Key exfilt_key_;

  // sub PKs
  Teuchos::RCP<ATS_Physics::PK_PhysicalBDF_Default> domain_flow_pk_;
  Teuchos::RCP<ATS_Physics::PK_PhysicalBDF_Default> surf_flow_pk_;

  // sub meshes
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;

  // coupled preconditioner
  Teuchos::RCP<Operators::Operator> precon_;
  Teuchos::RCP<Operators::Operator> precon_surf_;

  // Water delegate
  Teuchos::RCP<MPCDelegateWater> water_;
  bool consistent_cells_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> domain_db_;
  Teuchos::RCP<Debugger> surf_db_;

 private:
  // factory registration
  static RegisteredPKFactory<MPCCoupledWater> reg_;
};

} // namespace Amanzi


#endif
