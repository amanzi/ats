#include "heat_flux_bc_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

// registry of method
Utils::RegisteredFactory<FieldEvaluator,HeatFluxBCEvaluator> HeatFluxBCEvaluator::factory_("lake heat flux bc");

} //namespace
} //namespace
