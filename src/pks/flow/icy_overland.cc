/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "icy_overland.hh"

namespace Amanzi {
namespace Flow {

void IcyOverlandFlow::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // ensure that the overland conductivity uses the unfrozen ponded depth
  // -- set the height key to be eta * h, not just h, for the frozen case.
  //AMANZI_ASSERT(plist_->isSublist("overland conductivity evaluator") || plist_->isSublist("overland conductivity subgrid evaluator"));

  if (plist_->isSublist("overland conductivity evaluator")) {
    if (!plist_->sublist("overland conductivity evaluator").isParameter("depth key")) {
      plist_->sublist("overland conductivity evaluator").set("depth key",
              Keys::getKey(domain_,"unfrozen_effective_depth"));
    }
    AMANZI_ASSERT(plist_->sublist("overland conductivity evaluator").get<std::string>("depth key") != Keys::getKey(domain_,"ponded_depth"));
  }
  else if (plist_->isSublist("overland conductivity subgrid evaluator")) {
    if (!plist_->sublist("overland conductivity subgrid evaluator").isParameter("depth key")) {
      plist_->sublist("overland conductivity subgrid evaluator").set("depth key",
              Keys::getKey(domain_,"unfrozen_effective_depth"));
    }
    AMANZI_ASSERT(plist_->sublist("overland conductivity subgrid evaluator").get<std::string>("depth key") != Keys::getKey(domain_,"ponded_depth"));
  } 
  
  // Now continue as usual for overland head
  OverlandPressureFlow::SetupPhysicalEvaluators_(S);
}

} // namespace
} // namespace
