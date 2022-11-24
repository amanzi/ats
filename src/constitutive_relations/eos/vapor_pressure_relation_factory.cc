/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for Vapor Pressure implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "vapor_pressure_relation_factory.hh"

namespace Amanzi {
namespace Relations {

// method for instantiating VaporPressure implementations
Teuchos::RCP<VaporPressureRelation>
VaporPressureRelationFactory::createVaporPressure(Teuchos::ParameterList& plist)
{
  std::string eos_typename = plist.get<std::string>("vapor pressure model type");
  return Teuchos::rcp(CreateInstance(eos_typename, plist));
};

} // namespace Relations
} // namespace Amanzi
