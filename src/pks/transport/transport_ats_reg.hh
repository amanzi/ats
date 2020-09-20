/*
  Transport PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatskiy (dasvyat@lanl.gov)
*/

#include "transport_ats.hh"

namespace Amanzi {
namespace Transport {

RegisteredPKFactory<Transport_ATS> Transport_ATS::reg_("transport ATS");

}  // namespace Transport
}  // namespace Amanzi
