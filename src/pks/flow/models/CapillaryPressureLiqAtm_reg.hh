#include "CapillaryPressureLiqAtm.hh"

template<>
Amanzi::Utils::RegisteredFactory<Amanzi::Evaluator, Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::CapillaryPressureLiqAtm>>
Amanzi::EvaluatorModel_CompositeVector<ATS::Flow::Relations::CapillaryPressureLiqAtm>::Model_type::reg_("capillary pressure, atmospheric gas over liquid");

