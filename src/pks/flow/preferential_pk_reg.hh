/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
License: BSD
Authors: Daniil Svyatsky (dasvyat@lanl.gov)
         Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "preferential.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Preferential> Preferential::reg_("preferential flow");

} // namespace Flow
} // namespace Amanzi
