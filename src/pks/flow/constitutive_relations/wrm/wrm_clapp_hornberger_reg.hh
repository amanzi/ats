/*
  This is the flow component of the ATS physics code.
  License: BSD
  Authors: Markus Berndt (berndt@lanl.gov) 
  Konstantin Lipnikov (lipnikov@lanl.gov)
  F.-M. Yuan (yuanf@ornl.gov)
*/


#include "wrm_clapp_hornberger.hh"

namespace Amanzi {
namespace Flow {

Utils::RegisteredFactory<WRM,WRMClappHornberger> WRMClappHornberger::factory_("Clapp Hornberger");

}  // namespace
}  // namespace
