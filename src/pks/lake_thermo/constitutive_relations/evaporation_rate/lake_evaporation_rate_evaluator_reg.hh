#include "lake_evaporation_rate_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<Evaluator,LakeEvaporationRateEvaluator> LakeEvaporationRateEvaluator::factory_("lake evaporation rate");

} //namespace
} //namespace