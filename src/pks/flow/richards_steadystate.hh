/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Two-phase, variable density Richards equation in steady-state.
/*!

This is the same as Richards equation, but turns off the accumulation term.

.. _richards_steadystate-spec:
.. admonition:: richards_steadystate-spec

    INCLUDES:

    - ``[richards-spec]`` See `Richards PK`_

*/

#ifndef PK_FLOW_RICHARDS_STEADYSTATE_HH_
#define PK_FLOW_RICHARDS_STEADYSTATE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "State.hh"
#include "upwinding.hh"
#include "BoundaryFunction.hh"

#include "PK_Factory.hh"
#include "richards.hh"

namespace Amanzi {
namespace Flow {

class RichardsSteadyState : public Richards {
 public:
  // Constructors.

  RichardsSteadyState(Teuchos::ParameterList& FElist,
                      const Teuchos::RCP<Teuchos::ParameterList>& plist,
                      const Teuchos::RCP<State>& S,
                      const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~RichardsSteadyState() override {}

 protected:
  // ConstantTemperature is a BDFFnBase
  // computes the non-linear functional g = g(t,u,udot)
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g) override;

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

 private:
  // factory registration
  static RegisteredPKFactory<RichardsSteadyState> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
