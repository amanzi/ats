/*
  Transport PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  Self-registering factory for mulscale porosity models.
*/

#ifndef MULTISCALE_TRANSPORT_POROSITY_FACTORY_HH_
#define MULTISCALE_TRANSPORT_POROSITY_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "MultiscaleTransportPorosity.hh"

namespace Amanzi {
namespace Transport {

class MultiscaleTransportPorosityFactory : public Utils::Factory<MultiscaleTransportPorosity> {
 public:
  Teuchos::RCP<MultiscaleTransportPorosity> Create(Teuchos::ParameterList& plist);
};

}  // namespace Transport
}  // namespace Amanzi

#endif
