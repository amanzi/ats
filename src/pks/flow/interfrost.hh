/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*

A three-phase flow PK based on the Interfrost code comparison.

*/

#ifndef PK_FLOW_INTERFROST_HH_
#define PK_FLOW_INTERFROST_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "State.hh"
#include "upwinding.hh"
#include "BoundaryFunction.hh"

#include "PK.hh"
#include "PK_Factory.hh"

#include "permafrost.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

class Interfrost : public Permafrost {
 public:
  // Constructors.
  Interfrost(Teuchos::ParameterList& FElist,
             const Teuchos::RCP<Teuchos::ParameterList>& plist,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& solution)
    : PK(FElist, plist, S, solution), Permafrost(FElist, plist, S, solution)
  {}

  // Virtual destructor
  virtual ~Interfrost() override {}

  // -- accumulation term
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

 protected:
  virtual void AddAccumulation_(const Teuchos::Ptr<CompositeVector>& g) override;

  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_() override;

 private:
  // factory registration
  static RegisteredPKFactory<Interfrost> reg_;
};

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi

#endif
