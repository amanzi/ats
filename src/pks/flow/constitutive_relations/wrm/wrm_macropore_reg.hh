/*
  Copyright 2010-201x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

/*
  This is the flow component of the Amanzi code.
  License: BSD
  Authors: Scott Painter 
  Daniil Svyatsky
*/


#include "wrm_macropore.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM, WRMMacropore> WRMMacropore::factory_("macropore");

} // namespace Flow
} // namespace Amanzi
