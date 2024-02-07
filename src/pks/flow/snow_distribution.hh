/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Preferential distribution of snow precip in low-lying areas.
/*!

This PK is a heuristic PK that distributes incoming snow precipitation using a
diffusion wave equation.  Think of it as an analogue to overland flow -- it
effectively ensures that new snow "flows downhill," due to a uniformly random
direction and strength wind, and lands on the lowest lying areas.

Tweaking the snow-manning_coefficient lets you play with how uniform the snow
layer ends up.  Most of the parameters are set by your snow precipitation input
data interval.  The details of this are a bit tricky mathematically, and it may
take some fiddling with parameters to do this correctly if your data is not
daily (which all defaults are set for).

.. _snow-distribution-spec:
.. admonition:: snow-distribution-spec

    * `"distribution time`" ``[double]`` **86400.** Interval of snow precip input dataset. `[s]`
    * `"precipitation function`" ``[function-spec]`` Snow precipitation Function_ spec.

    * `"diffusion`" ``[pde-diffusion-spec]`` Diffusion drives the distribution.
      Typically we use finite volume here.  See PDE_Diffusion_

    * `"diffusion preconditioner`" ``[pde-diffusion-spec]`` Inverse of the
      above.  Likely only Jacobian term options are needed here, as the others
      default to the same as the `"diffusion`" list.  See PDE_Diffusion_.

    * `"inverse`" ``[inverse-typed-spec]`` Inverse_ method for the solve.

    Not typically provided by the user, defaults are good:

    * `"accumulation preconditioner`" ``[pde-accumulation-spec]`` See PDE_Accumulation_.


.. todo::
    For this PK, all variable root names are hard-coded.  This should get changed.

*/


#ifndef PK_FLOW_SNOW_DISTRIBUTION_HH_
#define PK_FLOW_SNOW_DISTRIBUTION_HH_

#include "upwinding.hh"

#include "Function.hh"

#include "Operator.hh"
#include "PDE_Diffusion.hh"
#include "PDE_Accumulation.hh"

#include "PK_Factory.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {
namespace Flow {

class SnowDistribution : public PK_PhysicalBDF_Default {
 public:
  SnowDistribution(Teuchos::ParameterList& pk_tree,
                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~SnowDistribution() {}

  // main methods
  // -- Initialize owned (dependent) variables.
  void Setup() override;

  // -- Initialize owned (dependent) variables.
  void Initialize() override;

  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  void FunctionalResidual(double t_old,
                          double t_new,
                          Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> g) override;

  // applies preconditioner to u and returns the result in Pu
  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  // updates the preconditioner
  void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // error monitor
  double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;

  bool
  ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override;

  // Choose a time step compatible with physics.
  double get_dt() override { return dt_factor_; }

  // Advance PK from time t_old to time t_new. True value of the last
  // parameter indicates drastic change of boundary and/or source terms
  // that may need PK's attention.
  //
  //  ALL SORTS OF FRAGILITY and UGLINESS HERE!
  //  DO NOT USE THIS OUT IN THE WILD!
  //
  //  1. this MUST go first
  //  2. it must be PERFECT NON_OVERLAPPING with everything else.  I'm not
  //     sure exactly what that means.  Something like, nothing that this PK
  //     writes can be read by anything else, except for the precip snow at
  //     the end?  But it should be ok?
  //  3. Exatrapolating in the timestepper should break things, so don't.
  //  4. set: pk's distribution time, potential's dt factor
  bool AdvanceStep(double t_old, double t_new, bool reinit) override;

  // -- Commit any secondary (dependent) variables.
  void CommitStep(double t_old, double t_new, const Tag& tag) override
  {
    // here to keep the coordinator from calling CommitSolution() since our
    // Advance() does it already
  }


 protected:
  // setup methods
  void SetupSnowDistribution_();
  void SetupPhysicalEvaluators_();

  // computational concerns in managing abs, rel perm
  // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  bool UpdatePermeabilityData_(const Tag& tag);

  // boundary condition members
  void UpdateBoundaryConditions_(const Tag& tag);

  // physical methods
  // -- diffusion term
  void ApplyDiffusion_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& g);
  // -- accumulation term
  void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g);


  Key precip_func_key_;
  Key cv_key_;
  Key cond_key_;
  Key elev_key_;
  Key precip_key_;
  Key uw_cond_key_;
  Key flux_dir_key_;
  Key potential_key_;
  // control switches
  Operators::UpwindMethod upwind_method_;

  double dt_factor_;
  double my_next_time_;

  // function for precip
  Teuchos::RCP<Function> precip_func_;

  // work data space
  Teuchos::RCP<Operators::Upwinding> upwinding_;

  // mathematical operators
  Teuchos::RCP<Operators::Operator> matrix_; // pc in PKPhysicalBDFBase
  Teuchos::RCP<Operators::PDE_Diffusion> matrix_diff_;
  Teuchos::RCP<Operators::PDE_Diffusion> face_matrix_diff_;
  Teuchos::RCP<Operators::PDE_Diffusion> preconditioner_diff_;
  Teuchos::RCP<Operators::PDE_Accumulation> preconditioner_acc_;

  // factory registration
  static RegisteredPKFactory<SnowDistribution> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
