#include "ManningConductivity.hh"

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::ManningConductivity>>
Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::ManningConductivity>::Model_type::reg_("Manning conductivity");


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

#include "CompressiblePorosityLinear.hh"

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModelByMaterial<ATS::Flow::Relations::CompressiblePorosityLinear>>
Amanzi::EvaluatorModelByMaterial<ATS::Flow::Relations::CompressiblePorosityLinear>::Model_type::by_material_reg_("compressible porosity by material");

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::CompressiblePorosityLinear>>
Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::CompressiblePorosityLinear>::Model_type::global_reg_("compressible porosity");

#include "CapillaryPressureLiqAtm.hh"

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::CapillaryPressureLiqAtm>>
Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::CapillaryPressureLiqAtm>::Model_type::reg_("capillary pressure, atmospheric gas over liquid");

#include "SurfaceWaterDepth.hh"

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::SurfaceWaterDepth>>
Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::SurfaceWaterDepth>::Model_type::reg_("surface water depth");

#include "SurfaceWaterContent.hh"
template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::SurfaceWaterContent>>
Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::SurfaceWaterContent>::Model_type::reg_("surface water content");

