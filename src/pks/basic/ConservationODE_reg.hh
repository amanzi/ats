#include "ConservationODE.hh"

template<>
Amanzi::RegisteredPKFactory<ATS::Basic::PK_ConservationODE_Implicit> ATS::Basic::PK_ConservationODE_Implicit::reg_("conservation ODE");

template<>
Amanzi::RegisteredPKFactory<ATS::Basic::PK_ConservationODE_Explicit> ATS::Basic::PK_ConservationODE_Explicit::reg_("conservation ODE, explicit");

template<>
Amanzi::RegisteredPKFactory<ATS::Basic::PK_ConservationODE_PredictorCorrector> ATS::Basic::PK_ConservationODE_PredictorCorrector::reg_("conservation ODE, predictor-corrector");
