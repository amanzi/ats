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

#include "transport_ats.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Transport {

RegisteredPKFactory<Transport_ATS> Transport_ATS::reg_("transport ATS");

} // namespace Transport
} // namespace ATS_Physics
} // namespace Amanzi
