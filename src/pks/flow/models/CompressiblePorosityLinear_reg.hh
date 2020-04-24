#include "CompressiblePorosityLinear.hh"

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModelByMaterial<ATS::Flow::Relations::CompressiblePorosityLinear>>
Amanzi::EvaluatorModelByMaterial<ATS::Flow::Relations::CompressiblePorosityLinear>::Model_type::by_material_reg_("compressible porosity by material");

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::CompressiblePorosityLinear>>
Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::CompressiblePorosityLinear>::Model_type::global_reg_("compressible porosity");

