/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)
*/

#include "preferential.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

RegisteredPKFactory<Preferential> Preferential::reg_("preferential flow");

} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
