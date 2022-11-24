/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#include <string>
#include "surface_relperm_model_factory.hh"

// explicity instantitate the static data of Factory<SurfaceRelPermModel>
template <>
Amanzi::Utils::Factory<Amanzi::Flow::SurfaceRelPermModel>::map_type*
  Amanzi::Utils::Factory<Amanzi::Flow::SurfaceRelPermModel>::map_;

namespace Amanzi {
namespace Flow {

// method for instantiating SurfaceRelPermModel implementations
Teuchos::RCP<SurfaceRelPermModel>
SurfaceRelPermModelFactory::createModel(Teuchos::ParameterList& plist)
{
  std::string type_name = plist.get<std::string>("surface rel perm model type");
  return Teuchos::rcp(CreateInstance(type_name, plist));
};

} // namespace Flow
} // namespace Amanzi
