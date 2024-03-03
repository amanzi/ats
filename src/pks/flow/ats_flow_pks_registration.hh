/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#pragma once

#include "richards.hh"
#include "richards_steadystate.hh"
#include "overland_pressure.hh"

namespace Amanzi {
namespace Flow {

REGISTER_PK(Richards);
REGISTER_PK(RichardsSteadyState);
REGISTER_PK(OverlandPressureFlow);

} // namespace Flow
} // namespace Amanzi
