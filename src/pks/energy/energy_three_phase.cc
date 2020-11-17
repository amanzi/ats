/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */
#include "primary_variable_field_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "thermal_conductivity_threephase_evaluator.hh"
#include "energy_three_phase.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void ThreePhase::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // Get data and evaluators needed by the PK
  // -- energy, the conserved quantity

  S->RequireField(conserved_key_)->SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator(conserved_key_);

  // -- thermal conductivity
  S->RequireField(conductivity_key_)->SetMesh(mesh_)
    ->SetGhosted()->AddComponent("cell", AmanziMesh::CELL, 1);
  Teuchos::ParameterList tcm_plist =
    plist_->sublist("thermal conductivity evaluator");
  tcm_plist.set("evaluator name", conductivity_key_);
  Teuchos::RCP<Energy::ThermalConductivityThreePhaseEvaluator> tcm =
    Teuchos::rcp(new Energy::ThermalConductivityThreePhaseEvaluator(tcm_plist));
  S->SetFieldEvaluator(conductivity_key_, tcm);

}

void
ThreePhase::Initialize(const Teuchos::Ptr<State>& S) {
  // INTERFROST comparison needs some very specialized ICs
  Teuchos::ParameterList& ic_plist = plist_->sublist("initial condition");
  if (ic_plist.isParameter("interfrost initial condition")) {
    std::string interfrost_ic = ic_plist.get<std::string>("interfrost initial condition");
    AMANZI_ASSERT(interfrost_ic == "TH3");

    Teuchos::RCP<CompositeVector> temp = S->GetFieldData(key_, name_);
    double r_sq = std::pow(0.5099,2);
    Epetra_MultiVector& temp_c = *temp->ViewComponent("cell", false);
    for (int c = 0; c!=temp_c.MyLength(); ++c) {
      AmanziGeometry::Point centroid = mesh_->cell_centroid(c);
      double circle_y = centroid[1] >= 0.5 ? 1.1 : -0.1;

      double dist = std::pow(centroid[0] - 0.5, 2) + std::pow(centroid[1] - circle_y, 2);
      if (dist <= r_sq) {
        temp_c[0][c] = 273.15 - 5.;
      } else {
        temp_c[0][c] = 273.15 + 5.;
      }
    }

    Teuchos::RCP<Field> field = S->GetField(key_, name_);
    field->set_initialized();

    // additionally call Initalize() to get faces from cell values
    field->Initialize(ic_plist);
  }

  TwoPhase::Initialize(S);
}


} // namespace Energy
} // namespace Amanzi
