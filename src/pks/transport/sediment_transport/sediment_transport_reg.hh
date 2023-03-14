/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy (dasvyat@lanl.gov)
*/

/*
  Transport PK

*/

#include "sediment_transport_pk.hh"

namespace Amanzi {
namespace SedimentTransport {

RegisteredPKFactory<SedimentTransport_PK> SedimentTransport_PK::reg_("sediment transport");

} // namespace SedimentTransport
} // namespace Amanzi
