/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "iem_quadratic.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<IEM, IEMQuadratic> IEMQuadratic::factory_("quadratic");

} // namespace Energy
} // namespace Amanzi
