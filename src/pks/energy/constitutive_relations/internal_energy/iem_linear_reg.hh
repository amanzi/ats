/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

------------------------------------------------------------------------- */

#include "iem_linear.hh"

namespace Amanzi {
namespace Energy {

Utils::RegisteredFactory<IEM, IEMLinear> IEMLinear::factory_("linear");

} // namespace Energy
} // namespace Amanzi
