/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! A three-phase, thermal Richard's equation with water, water vapor, and ice for permafrost applications.

/*!

Note that the only difference between permafrost and richards is in
constitutive relations -- the WRM changes to provide three saturations,
while the water content changes to account for water in ice phase.  As these
are now drop-in field evaluators, there is very little to change in the PK.

In the future, this should not even need a different PK.

.. _permafrost-spec:
.. admonition:: permafrost-spec
  
    * `"saturation ice key`" ``[string]`` **"DOMAIN-saturation_ice"** volume fraction of the ice phase (only when relevant) ``[-]`` Typically the default is correct.

    INCLUDES:

    - ``[richards-spec]`` See `Richards PK`_
*/

#ifndef PK_FLOW_PERMAFROST_HH_
#define PK_FLOW_PERMAFROST_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"
#include "TreeVector.hh"
#include "State.hh"
#include "upwinding.hh"
#include "BoundaryFunction.hh"

// #include "pk.hh"
// #include "pk_factory_ats.hh"
#include "PK.hh"
#include "PK_Factory.hh"

#include "richards.hh"

namespace Amanzi {
namespace Flow {

class Permafrost : public Richards {

friend class MPCCoupledFlowEnergy;

public:
  // Constructors.

  Permafrost(Teuchos::ParameterList& FElist,
             const Teuchos::RCP<Teuchos::ParameterList>& plist,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S, solution),
    Richards(FElist, plist, S, solution) {}

  // Virtual destructor
  virtual ~Permafrost() {}

protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S);

private:
  // factory registration
  static RegisteredPKFactory<Permafrost> reg_;

};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
