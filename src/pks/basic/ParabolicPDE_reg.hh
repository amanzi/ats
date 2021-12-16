#include "ParabolicPDE_MixedFormImplicit.hh"
template<>
Amanzi::RegisteredPKFactory<ATS::Basic::PK_ParabolicPDE_MixedFormImplicit> ATS::Basic::PK_ParabolicPDE_MixedFormImplicit::reg_("mixed form parabolic PDE");

#include "ParabolicPDE_PrimaryFormExplicit.hh"
template<>
Amanzi::RegisteredPKFactory<ATS::Basic::PK_ParabolicPDE_PrimaryFormExplicit> ATS::Basic::PK_ParabolicPDE_PrimaryFormExplicit::reg_("primary form parabolic PDE, explicit");

#include "ParabolicPDE_PredictorCorrector.hh"
template<>
Amanzi::RegisteredPKFactory<ATS::Basic::PK_ParabolicPDE_PredictorCorrector> ATS::Basic::PK_ParabolicPDE_PredictorCorrector::reg_("mixed form parabolic PDE, predictor-corrector");
