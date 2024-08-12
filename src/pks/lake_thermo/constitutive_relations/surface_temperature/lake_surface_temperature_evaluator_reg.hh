#include "lake_surface_temperature_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,LakeSurfaceTemperatureEvaluator> LakeSurfaceTemperatureEvaluator::factory_("lake surface temperature");

} //namespace
} //namespace
