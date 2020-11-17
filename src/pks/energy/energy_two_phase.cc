/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Process kernel for energy equation for Richard's flow.
------------------------------------------------------------------------- */


#include "eos_evaluator_tp.hh"
#include "iem_evaluator.hh"
#include "thermal_conductivity_twophase_evaluator.hh"
#include "enthalpy_evaluator.hh"
#include "energy_bc_factory.hh"

#include "energy_two_phase.hh"

namespace Amanzi {
namespace Energy {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------

TwoPhase::TwoPhase(Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution) :
    PK(FElist, plist, S, solution),
    EnergyBase(FElist, plist, S, solution) {}


// -------------------------------------------------------------
// Create the physical evaluators for energy, enthalpy, thermal
// conductivity, and any sources.
// -------------------------------------------------------------
void TwoPhase::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
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
  Teuchos::RCP<Energy::ThermalConductivityTwoPhaseEvaluator> tcm =
    Teuchos::rcp(new Energy::ThermalConductivityTwoPhaseEvaluator(tcm_plist));
  S->SetFieldEvaluator(conductivity_key_, tcm);

}


// -------------------------------------------------------------
// Initialize the needed models to plug in enthalpy.
// -------------------------------------------------------------
void TwoPhase::Initialize(const Teuchos::Ptr<State>& S) {
  // Call the base class's initialize.
  EnergyBase::Initialize(S);

  // For the boundary conditions, we currently hack in the enthalpy to
  // the boundary faces to correctly advect in a Dirichlet temperature
  // BC.  This requires density and internal energy, which in turn
  // require a model based on p,T.
  // This will be removed once boundary faces are implemented.
  Teuchos::RCP<FieldEvaluator> eos_fe =
    S->GetFieldEvaluator(Keys::getKey(domain_, "molar_density_liquid"));
  Teuchos::RCP<Relations::EOSEvaluator> eos_eval =
    Teuchos::rcp_dynamic_cast<Relations::EOSEvaluator>(eos_fe);
  AMANZI_ASSERT(eos_eval != Teuchos::null);
  eos_liquid_ = eos_eval->get_EOS();

  Teuchos::RCP<FieldEvaluator> iem_fe =
    S->GetFieldEvaluator(Keys::getKey(domain_, "internal_energy_liquid"));
  Teuchos::RCP<Energy::IEMEvaluator> iem_eval =
    Teuchos::rcp_dynamic_cast<Energy::IEMEvaluator>(iem_fe);
  AMANZI_ASSERT(iem_eval != Teuchos::null);
  iem_liquid_ = iem_eval->get_IEM();

}


// // -------------------------------------------------------------
// // Plug enthalpy into the boundary faces manually.
// // This will be removed once boundary faces exist.
// // -------------------------------------------------------------
// void TwoPhase::ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S) {
//   // put the boundary fluxes in faces for Dirichlet BCs.
//   // NOTE this boundary flux is in enthalpy, and
//   // h = n(T,p) * u_l(T) + p_l
//   Teuchos::RCP<const Epetra_MultiVector> pres;
//   Key pres_key = Keys::getKey(domain_, "pressure");
//   if (S->GetFieldData(pres_key)->HasComponent("face")) {
//     pres = S->GetFieldData(pres_key)->ViewComponent("face",false);
//   }
//   const Epetra_MultiVector& pres_c = *S->GetFieldData(pres_key)
//       ->ViewComponent("cell",false);

//   Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(key_);

//   const Epetra_MultiVector& flux = *S->GetFieldData(flux_key_)
//       ->ViewComponent("face",false);

//   bool include_work = plist_->sublist("enthalpy evaluator").get<bool>("include work term", true);
  
//   AmanziMesh::Entity_ID_List cells;
//   int nfaces = mesh_->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
//   for (int f=0; f!=nfaces; ++f) {
//     mesh_->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
//     if (bc_markers_adv_[f] == Operators::OPERATOR_BC_DIRICHLET) {
//       // If the advective markers are Dirichlet, and the diffusion markers are
//       // Neumann, that means we were given by the diffusive fluxes and the
//       // advected mass flux and temperature.
//       double T = bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET ? bc_values_[f] : BoundaryValue(temp, f);
//       double enthalpy = iem_liquid_->InternalEnergy(T);
//       if (include_work) {
//         double p = pres == Teuchos::null ? pres_c[0][cells[0]] : (*pres)[0][f];
//         double dens = eos_liquid_->MolarDensity(T,p);
//         enthalpy += p/dens;
//       }
//       bc_values_adv_[f] = enthalpy;
//     }
//   }
// }



} // namespace
} // namespace
