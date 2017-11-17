/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon


 <ParameterList name="boundary conditions">
   <ParameterList name="diffusive flux">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward diffusive flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


------------------------------------------------------------------------- */

#ifndef AMANZI_ENERGY_BC_FACTORY_HH_
#define AMANZI_ENERGY_BC_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "bc_factory.hh"

namespace Amanzi {
namespace Energy {

class EnergyBCFactory : public Amanzi::BCFactory {

public:
  EnergyBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                  Teuchos::ParameterList& plist) :
      Amanzi::BCFactory(mesh,plist) {}

  Teuchos::RCP<Functions::BoundaryFunction> CreateTemperature() const {
    return CreateWithFunction("temperature", "boundary temperature");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateTotalFlux() const {
    return CreateWithFunction("enthalpy flux", "outward enthalpy flux");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateDiffusiveFlux() const {
    return CreateWithFunction("diffusive flux", "outward diffusive flux");
  }
};

}  // namespace Energy
}  // namespace Amanzi

#endif
