/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------
ATS

Standard base for most diffusion-dominated PKs, this combines both
domains/meshes of PKPhysicalBase and Explicit methods of PKExplicitBase.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_EXPLICIT_DEFAULT_HH_
#define AMANZI_PK_PHYSICAL_EXPLICIT_DEFAULT_HH_

#include "errors.hh"
#include "PK.hh"
#include "pk_explicit_default.hh"
#include "PK_Physical.hh"


namespace Amanzi {


class PK_Physical_Explicit_Default : public PK_Explicit_Default, public PK_Physical {
 public:
  PK_Physical_Explicit_Default(Teuchos::ParameterList& pk_tree,
                               const Teuchos::RCP<Teuchos::ParameterList>& glist,
                               const Teuchos::RCP<State>& S,
                               const Teuchos::RCP<TreeVector>& solution)
    : PK(pk_tree, glist, S, solution),
      PK_Explicit_Default(pk_tree, glist, S, solution),
      PK_Physical(pk_tree, glist, S, solution)
  {}

  virtual void Setup() override
  {
    PK_Physical::Setup();
    PK_Explicit_Default::Setup();
  }

  // initialize.  Note both ExplicitBase and PhysicalBase have initialize()
  // methods, so we need a unique overrider.
  virtual void Initialize() override
  {
    PK_Physical::Initialize();
    PK_Explicit_Default::Initialize();
  }

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override
  {
    PK_Explicit_Default::AdvanceStep(t_old, t_new, reinit);
    ChangedSolutionPK(tag_next_);
    return false;
  }

  virtual void FailStep(double t_old, double t_new, const Tag& tag) override
  {
    PK_Physical::FailStep(t_old, t_new, tag);
  }
};

} // namespace Amanzi

#endif
