#include "MixedFormParabolicPDE.hh"

template<>
Amanzi::RegisteredPKFactory<ATS::Basic::PK_MixedFormParabolicPDE_Implicit> ATS::Basic::PK_MixedFormParabolicPDE_Implicit::reg_("mixed form parabolic PDE");

// template<>
// Amanzi::RegisteredPKFactory<ATS::Basic::PK_MixedFormParabolicPDE_Explicit> ATS::Basic::PK_MixedFormParabolicPDE_Explicit::reg_("mixed form parabolic PDE, explicit");

// template<>
// Amanzi::RegisteredPKFactory<ATS::Basic::PK_MixedFormParabolicPDE_PredictorCorrector> ATS::Basic::PK_MixedFormParabolicPDE_PredictorCorrector::reg_("mixed form parabolic PDE, predictor-corrector");
