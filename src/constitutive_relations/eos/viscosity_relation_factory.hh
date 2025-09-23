/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/* -------------------------------------------------------------------------

   ATS

   Self-registering factory for Viscosity implementations.
   ------------------------------------------------------------------------- */

#ifndef AMANZI_RELATIONS_VISCOSITY_RELATION_FACTORY_HH_
#define AMANZI_RELATIONS_VISCOSITY_RELATION_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "viscosity_relation.hh"
#include "Factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Relations {

class ViscosityRelationFactory : public Utils::Factory<ViscosityRelation> {
 public:
  Teuchos::RCP<ViscosityRelation> createViscosity(Teuchos::ParameterList& plist);
};

} // namespace Relations
} // namespace ATS_Physics
} // namespace Amanzi

#endif
