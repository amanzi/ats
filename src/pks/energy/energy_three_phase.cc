/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "EvaluatorPrimary.hh"
#include "enthalpy_evaluator.hh"
#include "thermal_conductivity_threephase_evaluator.hh"
#include "energy_three_phase.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
ThreePhase::ThreePhase(Teuchos::ParameterList& FElist,
                       const Teuchos::RCP<Teuchos::ParameterList>& plist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& solution)
  : PK(FElist, plist, S, solution), TwoPhase(FElist, plist, S, solution)
{}

// -------------------------------------------------------------
// Create the physical evaluators for energy and thermal conductivity
// -------------------------------------------------------------
void
ThreePhase::SetupPhysicalEvaluators_()
{
  // -- thermal conductivity
  // move evaluator from PK plist to State
  if (plist_->isSublist("thermal conductivity evaluator")) {
    auto& tcm_plist = S_->GetEvaluatorList(conductivity_key_);
    tcm_plist.setParameters(plist_->sublist("thermal conductivity evaluator"));
    tcm_plist.set("evaluator type", "three-phase thermal conductivity");
  }
  EnergyBase::SetupPhysicalEvaluators_();
}

void
ThreePhase::Initialize()
{
  // INTERFROST comparison needs some very specialized ICs
  if (plist_->sublist("initial condition").isParameter("interfrost initial condition")) {
    AMANZI_ASSERT(
      plist_->sublist("initial condition").get<std::string>("interfrost initial condition") ==
      "TH3");
    Teuchos::RCP<CompositeVector> temp = S_->GetPtrW<CompositeVector>(key_, tag_next_, name_);
    double r_sq = std::pow(0.5099, 2);
    Epetra_MultiVector& temp_c = *temp->ViewComponent("cell", false);
    for (int c = 0; c != temp_c.MyLength(); ++c) {
      AmanziGeometry::Point centroid = mesh_->getCellCentroid(c);
      double circle_y = centroid[1] >= 0.5 ? 1.1 : -0.1;

      double dist = std::pow(centroid[0] - 0.5, 2) + std::pow(centroid[1] - circle_y, 2);
      if (dist <= r_sq) {
        temp_c[0][c] = 273.15 - 5.;
      } else {
        temp_c[0][c] = 273.15 + 5.;
      }
    }

    S_->GetRecordW(key_, tag_next_, name_).set_initialized();

    // additionally call Initalize() to get faces from cell values
    S_->GetRecordW(key_, tag_next_, name_).Initialize(plist_->sublist("initial condition"));
  }

  EnergyBase::Initialize();
}


} // namespace Energy
} // namespace Amanzi
