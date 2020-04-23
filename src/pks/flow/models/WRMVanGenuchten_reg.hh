#include "WRMVanGenuchten.hh"

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModelByMaterial<ATS::Flow::Relations::WRMVanGenuchten>>
Amanzi::EvaluatorModelByMaterial<ATS::Flow::Relations::WRMVanGenuchten>::Model_type::by_material_reg_("WRM van Genuchten by material");

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::WRMVanGenuchten>>
Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::WRMVanGenuchten>::Model_type::global_reg_("WRM van Genuchten");

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModelByMaterial<ATS::Flow::Relations::WRMVanGenuchten_Kr>>
Amanzi::EvaluatorModelByMaterial<ATS::Flow::Relations::WRMVanGenuchten_Kr>::Model_type::by_material_reg_("Kr van Genuchten by material");

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::WRMVanGenuchten_Kr>>
Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::WRMVanGenuchten_Kr>::Model_type::global_reg_("Kr van Genuchten");

