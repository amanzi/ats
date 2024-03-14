/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_BC_FACTORY_HH_
#define AMANZI_FLOW_BC_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "bc_factory.hh"

/*!

Flow boundary conditions must follow the general format shown in
BoundaryConditions_.  Specific conditions implemented include:

Dirichlet (pressure) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for both surface and subsurface flows, this provides pressure data on
boundaries (in [Pa]).

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="pressure">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary pressure">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

Dirichlet (head) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Used for both surface and subsurface flows, this provides head data (in [m]
above the land surface), typically as a function of x & y.  In the subsurface
case, the z-value is given by hydrostatic relative to that head.

.. math::
  p = p_{atm} + rho * g * (h(x,y) + z_{surf} - z)

where h is the head function provided.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="head">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary head">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0.01"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Dirichlet (fixed level) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This fixes the water table at a constant elevation.  It is a head condition
that adapts to the surface elevation, adjusting the head to a datum that is a
fixed absolute z coordinate.

.. math::
  p = p_{atm} + rho * g * (h(x,y) - z)

where h is the head function provided.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="fixed level">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="fixed level">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Neumann (water flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for both surface and subsurface flows, this provides water flux data (in [mol m^-2 s^-1], for the subsurface, or [mol m^-1 s^-1] for the surface, in the outward normal direction) on boundaries.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="water flux">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward water flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="-1.e-3"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

Neumann (fix level flux) boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for surface only,  this provides fixed level ([m])  velocity data (in [m s^-1], in the outward normal direction on boundaries.

Example:

.. code-block:: xml

     <ParameterList name="boundary conditions">
       <ParameterList name="fixed level flux">
          <ParameterList name="river level south">
            <Parameter name="regions" type="Array(string)" value="{river south}"/>
            <ParameterList name="fixed level">
               <ParameterList name="function-constant">
                 <Parameter name="value" type="double" value="0.5"/>
               </ParameterList>
            </ParameterList>
            <ParameterList name="velocity">
               <ParameterList name="function-constant">
                 <Parameter name="value" type="double" value="2.5"/>
               </ParameterList>
            </ParameterList>
          </ParameterList>
       </ParameterList>
    </ParameterList>



Seepage face boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A variety of seepage face boundary conditions are permitted for both surface
and subsurface flow PKs.  Typically seepage conditions are of the form:

  - if :math:`q \cdot \hat{n} < 0`, then :math:`q = 0`
  - if :math:`p > p0`, then :math:`p = p0`

This ensures that flow is only out of the domain, but that the max pressure on
the boundary is specified by :math:`p0`.

Example: pressure (for surface or subsurface)

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="seepage face pressure">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary pressure">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Example: head (for surface)

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="seepage face head">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="boundary head">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Additionally, an infiltration flux may be prescribed, which describes the max
flux.  This is for surface faces on which a typical precipitation rate might
be prescribed, to be enforced until the water table rises to the surface, at
which point the precip is turned off and water seeps into runoff.  This
capability is experimental and has not been well tested.

  - if :math:`q \cdot \hat{n} < q_0`, then :math:`q = q_0`
  - if :math:`p > p_{atm}`, then :math:`p = p_{atm}`

Note the condition also accepts a parameter:

* `"explicit time index`" ``[bool]`` **false** If true, the _type_ of the BC is
  evaluated at the old time, keeping it fixed while the nonlinear solve
  iterates.

Example: seepage with infiltration

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="seepage face with infiltration">
     <ParameterList name="BC west">
       <Parameter name="explicit time index" type="bool" value="true"/>
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="outward water flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="-1.e-5"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
 </ParameterList>

Note it would be straightforward to add both p0 and q0 in the same condition;
this has simply not had a use case yet.


Zero head gradient boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Used for surface flows, this is an "outlet" boundary condition which looks to
enforce the condition that

.. math::
  \div h \cdot \hat{n} = 0

for head :math:`h` and outward normal :math:`\hat{n}`.  Note that this is an
"outlet" boundary, in the sense that it should really not be used on a
boundary in which

.. math::
  \div z \cdot \hat{n} > 0.

This makes it a useful boundary condition for benchmark and 2D problems, where
the elevation gradient is clear, but not so useful for DEM-based meshes.

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="zero gradient">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Critical depth boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Also for surface flows, this is an "outlet" boundary condition which looks to
set an outward flux to take away runoff.  This condition is given by:

.. math::
  q = \sqrt{g \hat{z}} n_{liq} h^1.5

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="critical depth">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Dynamic boundary condutions
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The type of boundary conditions maybe changed in time depending on the switch function of TIME.

.. code-block:: xml

   <ParameterList name="dynamic">
     <Parameter name="regions" type="Array(string)" value="{surface west}"/>
     <ParameterList name="switch function">
       <ParameterList name="function-tabular">
         <Parameter name="file" type="string" value="../data/floodplain2.h5" />
         <Parameter name="x header" type="string" value="Time" />
         <Parameter name="y header" type="string" value="Switch" />
         <Parameter name="form" type="Array(string)" value="{constant}"/>
       </ParameterList>
     </ParameterList>

     <ParameterList name="bcs">
       <Parameter name="bc types" type="Array(string)" value="{head, water flux}"/>
       <Parameter name="bc functions" type="Array(string)" value="{boundary head, outward water flux}"/>

       <ParameterList name="water flux">
         <ParameterList name="BC west">
           <Parameter name="regions" type="Array(string)" value="{surface west}"/>
           <ParameterList name="outward water flux">
             <ParameterList name="function-tabular">
               <Parameter name="file" type="string" value="../data/floodplain2.h5" />
               <Parameter name="x header" type="string" value="Time" />
               <Parameter name="y header" type="string" value="Flux" />
               <Parameter name="form" type="Array(string)" value="{linear}"/>
             </ParameterList>
            </ParameterList>
          </ParameterList>
       </ParameterList>

       <ParameterList name="head">
          <ParameterList name="BC west">
            <Parameter name="regions" type="Array(string)" value="{surface west}"/>
            <ParameterList name="boundary head">
              <ParameterList name="function-tabular">
                 <Parameter name="file" type="string" value="../data/floodplain2.h5" />
                 <Parameter name="x header" type="string" value="Time" />
                 <Parameter name="y header" type="string" value="Head" />
                 <Parameter name="form" type="Array(string)" value="{linear}"/>
               </ParameterList>
            </ParameterList>
          </ParameterList>
        </ParameterList>
     </ParameterList>

   </ParameterList>

 */


namespace Amanzi {
namespace Flow {

class FlowBCFactory : public Amanzi::BCFactory {
 public:
  FlowBCFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                const Teuchos::ParameterList& plist)
    : Amanzi::BCFactory(mesh, plist)
  {}

  Teuchos::RCP<Functions::BoundaryFunction> CreatePressure() const
  {
    return CreateWithFunction("pressure", "boundary pressure");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateHead() const
  {
    return CreateWithFunction("head", "boundary head");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateMassFlux() const
  {
    return CreateWithFunction("water flux", "outward water flux");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateZeroGradient() const
  {
    return CreateWithoutFunction("zero gradient");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateSeepageFaceHead() const
  {
    return CreateWithFunction("seepage face head", "boundary head");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateTidalHead() const
  {
    return CreateWithFunction("tidal head", "boundary head");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateSeepageFacePressure() const
  {
    return CreateWithFunction("seepage face pressure", "boundary pressure");
  }

  std::pair<bool, Teuchos::RCP<Functions::BoundaryFunction>>
  CreateSeepageFacePressureWithInfiltration()
  {
    bool is_explicit = CheckExplicitFlag("seepage face with infiltration");
    auto bfunc = CreateWithFunction("seepage face with infiltration", "outward water flux");
    return std::make_pair(is_explicit, bfunc);
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateCriticalDepth() const
  {
    return CreateWithoutFunction("critical depth");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateFixedLevel() const
  {
    return CreateWithFunction("fixed level", "fixed level");
  }

  Teuchos::RCP<Functions::DynamicBoundaryFunction> CreateDynamic() const
  {
    return CreateDynamicFunction("dynamic");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateFixedLevelFlux_Level() const
  {
    return CreateWithFunction("fixed level flux", "fixed level");
  }

  Teuchos::RCP<Functions::BoundaryFunction> CreateFixedLevelFlux_Velocity() const
  {
    return CreateWithFunction("fixed level flux", "velocity");
  }
};

} // namespace Flow
} // namespace Amanzi

#endif // AMANZI_FLOW_BC_FACTORY_HH_
