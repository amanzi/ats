/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  ATS

  EOS for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

*/

#include "errors.hh"
#include "viscosity_water.hh"

namespace Amanzi {
namespace Relations {

// registry of method
Utils::RegisteredFactory<ViscosityRelation, ViscosityWater> ViscosityWater::factory_(
  "liquid water");

} // namespace Relations
} // namespace Amanzi
