/* -*-  mode: c++; indent-tabs-mode: nil -*- */
// Boundary conditions base classes.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_BC_FACTORY_HH_
#define AMANZI_BC_FACTORY_HH_

/*!

In general, boundary conditions are provided in a hierarchical list by
boundary condition type, then functional form.  Boundary condition specs are
split between two types -- those which require a user-provided function
(i.e. Dirichlet data, etc) and those which do not (i.e. zero gradient
conditions).

A list of conditions might pull in both Dirichlet and Neumann data on
different regions, or use different functions on different regions.  The
following example illustrates how boundary conditions are prescribed across
the domain for a typical PK:

Example:

.. code-block:: xml

 <ParameterList name="boundary conditions">
   <ParameterList name="DIRICHLET_TYPE">
     <ParameterList name="BC west">
       <Parameter name="regions" type="Array(string)" value="{west}"/>
       <ParameterList name="DIRICHLET_FUNCTION_NAME">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="101325.0"/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
     <ParameterList name="BC east">
       <Parameter name="regions" type="Array(string)" value="{east}"/>
       <ParameterList name="DIRICHLET_FUNCTION_NAME">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="102325."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
   <ParameterList name="mass flux">
     <ParameterList name="BC north">
       <Parameter name="regions" type="Array(string)" value="{north}"/>
       <ParameterList name="outward mass flux">
         <ParameterList name="function-constant">
           <Parameter name="value" type="double" value="0."/>
         </ParameterList>
       </ParameterList>
     </ParameterList>
   </ParameterList>
   <ParameterList name="zero gradient">
     <ParameterList name="BC south">
       <Parameter name="regions" type="Array(string)" value="{south}"/>
     </ParameterList>
   </ParameterList>
 </ParameterList>


Different PKs populate this general format with different names, replacing
DIRICHLET_TYPE and DIRICHLET_FUNCTION_NAME.
  
 */


#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Point.hh"
#include "Mesh.hh"
#include "BoundaryFunction.hh"
#include "DynamicBoundaryFunction.hh"

namespace Amanzi {

class BCFactory {

public:
  BCFactory(const Teuchos::RCP<const AmanziMesh::Mesh> &mesh,
                const Teuchos::ParameterList& plist)
     : mesh_(mesh), plist_(plist) {}

  Teuchos::RCP<Functions::BoundaryFunction>
  CreateWithFunction(const std::string& list_name, const std::string& function_name) const;

  Teuchos::RCP<Functions::BoundaryFunction>
  CreateWithoutFunction(const std::string& list_name) const;

  Teuchos::RCP<Functions::DynamicBoundaryFunction>
  CreateDynamicFunction(const std::string& list_name) const;

  bool CheckExplicitFlag(const std::string& list_name);

 private:

  void ProcessListWithFunction_(const Teuchos::ParameterList&,
          const std::string& function_name,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessListWithoutFunction_(const Teuchos::ParameterList&,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessSpecWithFunction_(const Teuchos::ParameterList&,
          const std::string& function_name,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessSpecWithoutFunction_(const Teuchos::ParameterList&,
          const Teuchos::RCP<Functions::BoundaryFunction>&) const;

  void ProcessSpecWithFunctionRegions_(const Teuchos::ParameterList& list,
                                       const std::string& function_name,
                                       std::vector<std::string>& regions,
                                       const Teuchos::RCP<Functions::BoundaryFunction>& bc) const;

 private:
  const Teuchos::RCP<const AmanziMesh::Mesh>& mesh_;
  Teuchos::ParameterList plist_;
};

}  // namespace

#endif // AMANZI_BC_FACTORY_HH_
