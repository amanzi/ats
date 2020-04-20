#include "ConservationODE.hh"

template<>
Amanzi::RegisteredPKFactory<ATS::Basic::PK_ConservationODE_Implicit> ATS::Basic::PK_ConservationODE_Implicit::reg_("conservation ODE");
