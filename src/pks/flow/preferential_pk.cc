/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky(dasvyat@gmail.com)
*/

#include "boost/math/special_functions/fpclassify.hpp"

#include "boost/algorithm/string/predicate.hpp"

#include "Epetra_Import.h"

#include "flow_bc_factory.hh"

#include "Point.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "predictor_delegate_bc_flux.hh"
#include "wrm_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "richards_water_content_evaluator.hh"
#include "OperatorDefs.hh"
#include "BoundaryFlux.hh"
#include "pk_helpers.hh"

#include "preferential.hh"

#define DEBUG_RES_FLAG 0


namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------

Preferential::Preferential(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& glist,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution)
  : PK(pk_tree, glist, S, solution), Richards(pk_tree, glist, S, solution)
{
  coef_grav_key_ =
    Keys::readKey(*plist_, domain_, "gravity conductivity", "gravity_relative_permeability");

  // all manipulation of evaluator lists should happen in constructors (pre-setup)
  // -- Water retention evaluators for gravity term
  if (plist_->isSublist("water retention evaluator for gravity term")) {
    // note this overwrites the parameters from "water retention evaluator"
    // list set in Richards PK constructor
    auto& wrm_plist = S_->GetEvaluatorList(sat_key_);
    wrm_plist.setParameters(plist_->sublist("water retention evaluator for gravity term"));
  }
  if (S_->GetEvaluatorList(coef_grav_key_).numParams() == 0) {
    Teuchos::ParameterList& kr_plist = S_->GetEvaluatorList(coef_grav_key_);
    kr_plist.setParameters(S_->GetEvaluatorList(sat_key_));
    kr_plist.set<std::string>("evaluator type", "WRM rel perm");
  }
  S_->GetEvaluatorList(coef_grav_key_).set<double>("permeability rescaling", perm_scale_);
}


// -------------------------------------------------------------
// Pieces of the construction process that are common to all
// Preferential-like PKs.
// -------------------------------------------------------------
void
Preferential::RequireNonlinearCoefficient_(const Key& key, const std::string& coef_location)
{
  // -- require the data on appropriate locations
  if (coef_location == "upwind: face") {
    S_->Require<CompositeVector, CompositeVectorSpace>(key, tag_next_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponents({ "face", "grav" }, { AmanziMesh::Entity_kind::FACE, AmanziMesh::Entity_kind::FACE }, { 1, 1 });
  } else if (coef_location == "standard: cell") {
    S_->Require<CompositeVector, CompositeVectorSpace>(key, tag_next_, name_)
      .SetMesh(mesh_)
      ->SetGhosted()
      ->SetComponents({ "cell", "grav" }, { AmanziMesh::Entity_kind::CELL, AmanziMesh::Entity_kind::FACE }, { 1, 1 });
  } else {
    Errors::Message message("Unknown upwind coefficient location in Preferential flow.");
    Exceptions::amanzi_throw(message);
  }
  S_->GetRecordW(key, tag_next_, name_).set_io_vis(false);
}


// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void
Preferential::SetupPhysicalEvaluators_()
{
  Richards::SetupPhysicalEvaluators_();

  // -- rel perm for gravity term
  requireAtNext(coef_grav_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
}


// -----------------------------------------------------------------------------
// Use the physical rel perm (on cells) to update a work vector for rel perm.
//
//   This deals with upwinding, etc.
// -----------------------------------------------------------------------------
bool
Preferential::UpdatePermeabilityData_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating permeability?";

  Teuchos::RCP<const CompositeVector> rel_perm = S_->GetPtr<CompositeVector>(coef_key_, tag);
  Teuchos::RCP<const CompositeVector> rel_perm_grav =
    S_->GetPtr<CompositeVector>(coef_grav_key_, tag);
  bool update_perm = S_->GetEvaluator(coef_key_, tag).Update(*S_, name_);
  update_perm |= S_->GetEvaluator(coef_grav_key_, tag).Update(*S_, name_);

  // requirements due to the upwinding method
  if (Krel_method_ == Operators::UPWIND_METHOD_TOTAL_FLUX) {
    bool update_dir = S_->GetEvaluator(mass_dens_key_, tag).Update(*S_, name_);
    update_dir |= S_->GetEvaluator(key_, tag).Update(*S_, name_);

    if (update_dir) {
      // update the direction of the flux -- note this is NOT the flux
      Teuchos::RCP<const CompositeVector> rho = S_->GetPtr<CompositeVector>(mass_dens_key_, tag);
      Teuchos::RCP<CompositeVector> flux_dir =
        S_->GetPtrW<CompositeVector>(flux_dir_key_, tag, name_);
      Teuchos::RCP<const CompositeVector> pres = S_->GetPtr<CompositeVector>(key_, tag);

      if (!deform_key_.empty() &&
          S_->GetEvaluator(deform_key_, tag_next_).Update(*S_, name_ + " flux dir"))
        face_matrix_diff_->SetTensorCoefficient(K_);

      face_matrix_diff_->SetDensity(rho);
      face_matrix_diff_->UpdateMatrices(Teuchos::null, pres.ptr());
      //if (!pres->HasComponent("face"))
      face_matrix_diff_->ApplyBCs(true, true, true);
      face_matrix_diff_->UpdateFlux(pres.ptr(), flux_dir.ptr());

      if (clobber_boundary_flux_dir_) {
        Epetra_MultiVector& flux_dir_f = *flux_dir->ViewComponent("face", false);

        auto& markers = bc_markers();
        auto& values = bc_values();

        for (int f = 0; f != markers.size(); ++f) {
          if (markers[f] == Operators::OPERATOR_BC_NEUMANN) {
            auto cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
            AMANZI_ASSERT(cells.size() == 1);
            int c = cells[0];
            const auto& [faces, dirs] = mesh_->getCellFacesAndDirections(c);
            int i = std::find(faces.begin(), faces.end(), f) - faces.begin();

            flux_dir_f[0][f] = values[f] * dirs[i];
          }
        }
      }
    }

    update_perm |= update_dir;
  }

  if (update_perm) {
    Teuchos::RCP<CompositeVector> uw_rel_perm =
      S_->GetPtrW<CompositeVector>(uw_coef_key_, tag, name_);

    // // Move rel perm on boundary_faces into uw_rel_perm on faces
    const Epetra_Import& vandelay = mesh_->getBoundaryFaceImporter();
    const Epetra_MultiVector& rel_perm_bf = *rel_perm->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& rel_perm_grav_bf =
      *rel_perm_grav->ViewComponent("boundary_face", false);

    Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face", false);
    uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);
    Epetra_MultiVector& uw_rel_perm_grav = *uw_rel_perm->ViewComponent("grav", false);
    uw_rel_perm_grav.Export(rel_perm_grav_bf, vandelay, Insert);

    // Upwind, only overwriting boundary faces if the wind says to do so.
    upwinding_->Update(*rel_perm, "cell", *uw_rel_perm, "face", *S_);
    upwinding_->Update(*rel_perm, "cell", *uw_rel_perm, "grav", *S_);

    if (clobber_policy_ == "clobber") {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face", false);
      uw_rel_perm_f.Export(rel_perm_bf, vandelay, Insert);
    } else if (clobber_policy_ == "max") {
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face", false);
      const auto& fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
      const auto& bfmap = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
      for (int bf = 0; bf != rel_perm_bf.MyLength(); ++bf) {
        auto f = fmap.LID(bfmap.GID(bf));
        if (rel_perm_bf[0][bf] > uw_rel_perm_f[0][f]) { uw_rel_perm_f[0][f] = rel_perm_bf[0][bf]; }
      }
    } else if (clobber_policy_ == "unsaturated") {
      // clobber only when the interior cell is unsaturated
      Epetra_MultiVector& uw_rel_perm_f = *uw_rel_perm->ViewComponent("face", false);
      const Epetra_MultiVector& pres =
        *S_->Get<CompositeVector>(key_, tag).ViewComponent("cell", false);
      const auto& fmap = mesh_->getMap(AmanziMesh::Entity_kind::FACE,true);
      const auto& bfmap = mesh_->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,true);
      for (int bf = 0; bf != rel_perm_bf.MyLength(); ++bf) {
        auto f = fmap.LID(bfmap.GID(bf));
        auto fcells = mesh_->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        AMANZI_ASSERT(fcells.size() == 1);
        if (pres[0][fcells[0]] < 101225.) {
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf];
        } else if (pres[0][fcells[0]] < 101325.) {
          double frac = (101325. - pres[0][fcells[0]]) / 100.;
          uw_rel_perm_f[0][f] = rel_perm_bf[0][bf] * frac + uw_rel_perm_f[0][f] * (1 - frac);
        }
      }
    }

    if (uw_rel_perm->HasComponent("face")) uw_rel_perm->ScatterMasterToGhosted("face");
    if (uw_rel_perm->HasComponent("grav")) uw_rel_perm->ScatterMasterToGhosted("grav");
  }

  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) { *vo_->os() << " " << update_perm << std::endl; }

  return update_perm;
};


bool
Preferential::UpdatePermeabilityDerivativeData_(const Tag& tag)
{
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) *vo_->os() << "  Updating permeability derivatives?";

  bool update_perm = S_->GetEvaluator(coef_key_, tag).UpdateDerivative(*S_, name_, key_, tag);
  update_perm |= S_->GetEvaluator(coef_grav_key_, tag).UpdateDerivative(*S_, name_, key_, tag);
  ;

  if (update_perm) {
    const CompositeVector& drel_perm =
      S_->GetDerivative<CompositeVector>(coef_key_, tag, key_, tag);
    const CompositeVector& drel_grav_perm =
      S_->GetDerivative<CompositeVector>(coef_grav_key_, tag, key_, tag);

    if (!duw_coef_key_.empty()) {
      CompositeVector& duw_rel_perm = S_->GetW<CompositeVector>(duw_coef_key_, tag, name_);
      duw_rel_perm.PutScalar(0.);

      // Upwind, only overwriting boundary faces if the wind says to do so.
      upwinding_deriv_->Update(drel_perm, "cell", duw_rel_perm, "face", *S_);
      upwinding_deriv_->Update(drel_grav_perm, "cell", duw_rel_perm, "grav", *S_);

      duw_rel_perm.ScatterMasterToGhosted("face");
      duw_rel_perm.ScatterMasterToGhosted("grav");
    } else {
      drel_perm.ScatterMasterToGhosted("cell");
      drel_grav_perm.ScatterMasterToGhosted("cell");
    }
  }

  // debugging
  if (vo_->os_OK(Teuchos::VERB_EXTREME)) { *vo_->os() << " " << update_perm << std::endl; }
  return update_perm;
};

} // namespace Flow
} // namespace Amanzi
