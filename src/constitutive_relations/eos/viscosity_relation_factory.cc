/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include <string>
#include "viscosity_relation_factory.hh"

namespace Amanzi {
namespace Relations {

// method for instantiating Viscosity implementations
Teuchos::RCP<ViscosityRelation>
ViscosityRelationFactory::createViscosity(Teuchos::ParameterList& plist)
{
  std::string visc_typename = plist.get<std::string>("viscosity type");
  return Teuchos::rcp(CreateInstance(visc_typename, plist));
};

} // namespace Relations
} // namespace Amanzi
