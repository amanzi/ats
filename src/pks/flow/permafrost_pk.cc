/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "flow_bc_factory.hh"

#include "upwind_cell_centered.hh"
#include "upwind_arithmetic_mean.hh"
#include "upwind_total_flux.hh"
#include "upwind_gravity_flux.hh"

#include "CompositeVectorFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

#include "EvaluatorPrimary.hh"
#include "wrm_permafrost_evaluator.hh"
#include "rel_perm_evaluator.hh"
#include "rel_perm_sutraice_evaluator.hh"
#include "rel_perm_frzBC_evaluator.hh"
#include "pk_helpers.hh"

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void
Permafrost::SetupPhysicalEvaluators_()
{
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  requireAtNext(perm_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, num_perm_vals_);

  // -- water content, and evaluator, and derivative for PC
  requireAtNext(conserved_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    conserved_key_, tag_next_, key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  requireAtCurrent(conserved_key_, tag_current_, *S_, name_, true);

  // -- Water retention evaluators
  // This deals with deprecated location for the WRM list (in the PK).  Move it
  // to state.
  // -- This setup is a little funky -- we use four evaluators to capture the physics.
  if (plist_->isSublist("water retention evaluator")) {
    auto& wrm_plist = S_->GetEvaluatorList(sat_key_);
    wrm_plist.setParameters(plist_->sublist("water retention evaluator"));
    wrm_plist.set("evaluator type", "permafrost WRM");
  }
  if (!S_->HasEvaluator(coef_key_, tag_next_) &&
      (S_->GetEvaluatorList(coef_key_).numParams() == 0)) {
    Teuchos::ParameterList& kr_plist = S_->GetEvaluatorList(coef_key_);
    kr_plist.setParameters(S_->GetEvaluatorList(sat_key_));
    kr_plist.set<std::string>("evaluator type", "WRM rel perm");
  }

  // -- saturation
  requireAtNext(sat_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
  requireAtNext(sat_gas_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
  requireAtNext(sat_ice_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);
  auto& wrm = S_->RequireEvaluator(sat_key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  requireAtCurrent(sat_key_, tag_current_, *S_, name_, true);
  // S_->RequireEvaluator(sat_key_, tag_current_);
  requireAtCurrent(sat_ice_key_, tag_current_, *S_, name_, true);
  // S_->RequireEvaluator(sat_ice_key_, tag_current_);

  // -- the rel perm evaluator, also with the same underlying WRM.
  S_->GetEvaluatorList(coef_key_).set<double>("permeability rescaling", perm_scale_);
  requireAtNext(coef_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1)
    ->AddComponent("boundary_face", AmanziMesh::Entity_kind::BOUNDARY_FACE, 1);

  // -- get the WRM models
  auto wrm_eval = dynamic_cast<Flow::WRMPermafrostEvaluator*>(&wrm);
  AMANZI_ASSERT(wrm_eval != nullptr);
  wrms_ = wrm_eval->get_WRMs();

  // -- molar density used to infer liquid Darcy velocity from flux
  requireAtNext(molar_dens_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  // -- liquid mass density for the gravity fluxes
  requireAtNext(mass_dens_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
}

} // namespace Flow
} // namespace Amanzi
