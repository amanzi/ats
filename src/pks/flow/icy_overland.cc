/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */
#include "independent_variable_field_evaluator.hh"

#include "upwinding.hh"
#include "upwind_potential_difference.hh"
#include "pres_elev_evaluator.hh"
#include "elevation_evaluator.hh"
#include "meshed_elevation_evaluator.hh"
#include "standalone_elevation_evaluator.hh"
#include "overland_conductivity_evaluator.hh"
#include "overland_conductivity_model.hh"
#include "unfrozen_effective_depth_evaluator.hh"
#include "unfrozen_fraction_model.hh"
#include "unfrozen_fraction_evaluator.hh"

#include "overland_pressure_water_content_evaluator.hh"
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
