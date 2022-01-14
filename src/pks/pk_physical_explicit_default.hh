/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Standard base for most diffusion-dominated PKs, this combines both
domains/meshes of PKPhysicalBase and Explicit methods of PKExplicitBase.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_EXPLICIT_DEFAULT_HH_
#define AMANZI_PK_PHYSICAL_EXPLICIT_DEFAULT_HH_

#include "errors.hh"
#include "PK.hh"
#include "pk_explicit_default.hh"
#include "pk_physical_default.hh"


namespace Amanzi {


class  PK_Physical_Explicit_Default : public PK_Explicit_Default, public PK_Physical_Default {

public:
  PK_Physical_Explicit_Default(Teuchos::ParameterList& pk_tree,
                          const Teuchos::RCP<Teuchos::ParameterList>& glist,
                          const Teuchos::RCP<State>& S,
                          const Teuchos::RCP<TreeVector>& solution)
    : PK(pk_tree, glist, S, solution),
      PK_Explicit_Default(pk_tree, glist, S, solution),
      PK_Physical_Default(pk_tree, glist, S, solution) {}

  virtual void Setup(const Teuchos::Ptr<State>& S) {
    PK_Physical_Default::Setup(S);
    PK_Explicit_Default::Setup(S);
  }

  // initialize.  Note both ExplicitBase and PhysicalBase have initialize()
  // methods, so we need a unique overrider.
  virtual void Initialize(const Teuchos::Ptr<State>& S) {
    PK_Physical_Default::Initialize(S);
    PK_Explicit_Default::Initialize(S);
  }

  // -- Advance from state S0 to state S1 at time S0.time + dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) {
    PK_Explicit_Default::AdvanceStep(t_old, t_new, reinit);
    ChangedSolution();
    return false;
  }


  // -- Experimental approach -- calling this indicates that the time
  //    integration scheme is changing the value of the solution in
  //    state.
  virtual void ChangedSolution() {
    solution_evaluator_->SetChanged(S_next_.ptr());
  }

};

}

#endif
