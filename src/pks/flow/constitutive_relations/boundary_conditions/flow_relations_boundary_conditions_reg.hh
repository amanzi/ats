#include "Factory.hh"

#include "registration_macro.hh"

#include "EvaluatorBCPondedDepth.hh"
#include "EvaluatorBCSeepageFaceHead.hh"

namespace Amanzi {

REGISTER(Flow::Relations::EvaluatorBCPondedDepth);
REGISTER(Flow::Relations::EvaluatorBCSeepageFaceHead);


} // namespace Amanzi
