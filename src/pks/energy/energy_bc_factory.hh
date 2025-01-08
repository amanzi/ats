/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

Energy boundary conditions must follow the general format shown in
BoundaryConditions_.  Energy-specific conditions implemented include:

Dirichlet (temperature) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Provide temperature data in units of [K].

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="temperature">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary temperature">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="276.15."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Neumann (diffusive energy flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Note that for an energy equation, there are two mechanisms by which energy can
be fluxed into the domain: through diffusion and through advection.  This
boundary condition sets the diffusive flux of energy, and allows the advective
flux to be whatever it may be.  Frequently this is used in combination with
boundaries where water is expected to be advected out of the domain, and we
wish to allow the energy of that water to be advected away with it, but wish to
independently specify diffusive fluxes.  This can also be used in cases where
the mass flux is prescribed to be zero (e.g. bottom boundaries, where this
might be the geothermal gradient).

Units are in **[MW m^-2]**, noting the deviation from SI units!

Example:

.. code-block:: xml

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

Note that another commonly implemented boundary condition is one where the
diffusive flux is prescribed, and also the temperature of incoming water is
prescribed.  This is not currently implemented, but would be straightforward to
do so if requested.


Neumann (total energy flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This boundary condition sets the total flux of energy, from both advection and
diffusion.  This is not used all that often in real applications, but is common
for benchmarks or other testing.

Units are in **[MW m^-2]**, noting the deviation from SI units!

Example:

.. code-block:: xml

  <ParameterList name="boundary conditions">
   <ParameterList name="enthalpy flux">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward enthalpy flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>



*/

#ifndef AMANZI_ENERGY_BC_FACTORY_HH_
#define AMANZI_ENERGY_BC_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "bc_factory.hh"

namespace Amanzi {
namespace Energy {

class EnergyBCFactory : public Amanzi::BCFactory {
 public:
  EnergyBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh, Teuchos::ParameterList& plist)
    : Amanzi::BCFactory(mesh, plist)
  {}

  Teuchos::RCP<Functions::BoundaryFunction> CreateTemperature() const
  {
    return CreateWithFunction("temperature", "boundary temperature");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateTotalFlux() const
  {
    return CreateWithFunction("enthalpy flux", "outward enthalpy flux");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateDiffusiveFlux() const
  {
    return CreateWithFunction("diffusive flux", "outward diffusive flux");
  }
};

} // namespace Energy
} // namespace Amanzi

#endif
