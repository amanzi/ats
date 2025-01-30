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
  // -- water content, and evaluator, and derivative for PC
  requireAtNext(conserved_key_, tag_next_, *S_)
    .SetMesh(mesh_)
    ->SetGhosted()
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S_->RequireDerivative<CompositeVector, CompositeVectorSpace>(
    conserved_key_, tag_next_, key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  requireAtCurrent(conserved_key_, tag_current_, *S_, name_, true);

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
  requireAtCurrent(sat_ice_key_, tag_current_, *S_, name_, true);

  // -- the rel perm evaluator, also with the same underlying WRM.
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
