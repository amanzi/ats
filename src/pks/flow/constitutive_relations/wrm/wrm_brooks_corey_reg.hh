/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Bo Gao (gaob@ornl.gov)
*/

#include "wrm_brooks_corey.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMBrooksCorey> WRMBrooksCorey::factory_("Brooks-Corey");

} // namespace Flow
} // namespace Amanzi
