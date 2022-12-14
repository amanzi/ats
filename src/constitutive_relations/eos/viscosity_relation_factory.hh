/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for Viscosity implementations.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_RELATIONS_VISCOSITY_RELATION_FACTORY_HH_
#define AMANZI_RELATIONS_VISCOSITY_RELATION_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "viscosity_relation.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Relations {

class ViscosityRelationFactory : public Utils::Factory<ViscosityRelation> {
 public:
  Teuchos::RCP<ViscosityRelation> createViscosity(Teuchos::ParameterList& plist);
};

} // namespace Relations
} // namespace Amanzi

#endif
