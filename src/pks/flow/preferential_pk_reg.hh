/* -*-  mode: c++; indent-tabs-mode: nil -*- */

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

} // namespace
} // namespace
