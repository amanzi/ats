#include "EvaluatorModelCV.hh"
#include "EvaluatorModelCVByMaterial.hh"
#include "Factory.hh"

#include "registration_macro.hh"
#include "soil_resistance_sakagucki_zeng_evaluator.hh"

namespace Amanzi {

template <>
REGISTER(Flow::Relations::SoilResistanceSakaguckiZengEvaluator);
template <>
REGISTER(Flow::Relations::SoilResistanceSakaguckiZengEvaluatorByMaterial);

} // namespace Amanzi
