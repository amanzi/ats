#include "richards.hh"


Amanzi::RegisteredPKFactory<ATS::Flow::Richards_Implicit> ATS::Flow::Richards_Implicit::reg_("richards flow, implicit");
Amanzi::RegisteredPKFactory<ATS::Flow::Richards_Explicit> ATS::Flow::Richards_Explicit::reg_("richards flow, explicit");

