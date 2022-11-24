/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatskiy (dasvyat@lanl.gov)
*/

#include "sediment_transport_pk.hh"

namespace Amanzi {
namespace SedimentTransport {

RegisteredPKFactory<SedimentTransport_PK> SedimentTransport_PK::reg_("sediment transport");

} // namespace SedimentTransport
} // namespace Amanzi
