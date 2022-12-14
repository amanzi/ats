/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code.
------------------------------------------------------------------------- */

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<Permafrost> Permafrost::reg_("permafrost flow");

} // namespace Flow
} // namespace Amanzi
