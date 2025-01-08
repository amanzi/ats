/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "overland_pressure.hh"

namespace Amanzi {
namespace Flow {

RegisteredPKFactory<OverlandPressureFlow>
  OverlandPressureFlow::reg_("overland flow, pressure basis");


} // namespace Flow
} // namespace Amanzi
