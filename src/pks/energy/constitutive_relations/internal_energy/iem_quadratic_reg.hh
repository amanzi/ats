/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "iem_quadratic.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

Utils::RegisteredFactory<IEM, IEMQuadratic> IEMQuadratic::factory_("quadratic");

} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
