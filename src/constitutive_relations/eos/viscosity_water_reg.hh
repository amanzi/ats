/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "errors.hh"
#include "viscosity_water.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<ViscosityRelation, ViscosityWater>
  ViscosityWater::factory_("liquid water");

} // namespace Relations
} // namespace Amanzi
