/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

   ATS
   Author: Ethan Coon

   Self-registering factory for EOS implementations.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_RELATIONS_VAPOR_PRESSURE_MODEL_FACTORY_HH_
#define AMANZI_RELATIONS_VAPOR_PRESSURE_MODEL_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "vapor_pressure_relation.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Relations {

class VaporPressureRelationFactory : public Utils::Factory<VaporPressureRelation> {
 public:
  Teuchos::RCP<VaporPressureRelation> createVaporPressure(Teuchos::ParameterList& plist);
};

} // namespace Relations
} // namespace Amanzi

#endif
