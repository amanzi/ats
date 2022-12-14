/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "richards.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Richards> Richards::reg_("richards flow");

} // namespace Flow
} // namespace Amanzi
