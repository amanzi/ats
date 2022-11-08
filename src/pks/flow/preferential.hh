/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: 
*/
//! 

/*!

Solves Preferential Flow equation:


*/

#ifndef PK_FLOW_PREFERENTIAL_HH_
#define PK_FLOW_PREFERENTIAL_HH_

#include "richards.hh"

#include "PDE_DiffusionFactory.hh"
#include "PDE_Accumulation.hh"
#include "PK_Factory.hh"
#include "pk_physical_bdf_default.hh"

namespace Amanzi {

// forward declarations
class MPCSubsurface;
class PredictorDelegateBCFlux;
class PrimaryVariableFieldEvaluator;
namespace WhetStone { class Tensor; }

namespace Flow {

class Preferential : public Richards {

public:

  Preferential(Teuchos::ParameterList& pk_tree,
           const Teuchos::RCP<Teuchos::ParameterList>& plist,
           const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~Preferential() {}

  // -- Setup data.
  virtual void Setup() override;


protected:
  // Create of physical evaluators.
  virtual void SetupPhysicalEvaluators_() override;
  virtual void SetupPreferentialFlow_();

  // // boundary condition members
  // void ComputeBoundaryConditions_(const Teuchos::Ptr<State>& S);

  // // -- builds tensor K, along with faced-based Krel if needed by the rel-perm method
  virtual bool UpdatePermeabilityData_(const Tag& tag) override;
  virtual bool UpdatePermeabilityDerivativeData_(const Tag& tag) override;


protected:
   Key coef_grav_key_, dcoef_grav_key_; 
  
private:
  // factory registration
  static RegisteredPKFactory<Preferential> reg_;

  // Preferential has a friend in couplers...
  friend class Amanzi::MPCSubsurface;

};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
