/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   ATS

   Self-registering factory for Vapor Pressure implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "vapor_pressure_relation_factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

// method for instantiating VaporPressure implementations
Teuchos::RCP<VaporPressureRelation>
VaporPressureRelationFactory::createVaporPressure(Teuchos::ParameterList& plist)
{
  std::string eos_typename = plist.get<std::string>("vapor pressure model type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi
