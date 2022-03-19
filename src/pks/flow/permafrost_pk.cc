/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
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
#include "wrm_evaluator.hh"
#include "richards_water_content_evaluator.hh"
#include "pk_helpers.hh"

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Permafrost::SetupPhysicalEvaluators_()
{
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S_->Require<CompositeVector,CompositeVectorSpace>(perm_key_, tag_next_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, num_perm_vals_);
  S_->RequireEvaluator(perm_key_, tag_next_);

  // -- water content, and evaluator, and derivative for PC
  S_->Require<CompositeVector,CompositeVectorSpace>(conserved_key_, tag_next_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireDerivative<CompositeVector,CompositeVectorSpace>(conserved_key_,
          tag_next_, key_, tag_next_);
  S_->RequireEvaluator(conserved_key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  S_->Require<CompositeVector,CompositeVectorSpace>(conserved_key_, tag_current_, name_);
  // RequireEvaluatorPrimary(conserved_key_, tag_current_, *S_);

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
  S_->Require<CompositeVector,CompositeVectorSpace>(sat_key_, tag_next_).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->Require<CompositeVector,CompositeVectorSpace>(sat_gas_key_, tag_next_).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1) 
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->Require<CompositeVector,CompositeVectorSpace>(sat_ice_key_, tag_next_).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);

  // -- capillary pressure
  Key capillary_pressure_gas_liq_key_ = Keys::readKey(*plist_, domain_, 
      "capillary_pressure_gas_liq", "capillary_pressure_gas_liq");
  Key capillary_pressure_liq_ice_key_ = Keys::readKey(*plist_, domain_, 
      "capillary_pressure_liq_ice", "capillary_pressure_liq_ice");
  S_->Require<CompositeVector,CompositeVectorSpace>(capillary_pressure_gas_liq_key_, tag_next_)
      .SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->Require<CompositeVector,CompositeVectorSpace>(capillary_pressure_liq_ice_key_, tag_next_)
      .SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1)
      ->AddComponent("boundary_face", AmanziMesh::BOUNDARY_FACE, 1);
  S_->RequireEvaluator(capillary_pressure_gas_liq_key_, tag_next_);
  S_->RequireEvaluator(capillary_pressure_liq_ice_key_, tag_next_);

   auto& wrm = S_->RequireEvaluator(sat_key_, tag_next_);
  S_->RequireEvaluator(sat_gas_key_, tag_next_);
  S_->RequireEvaluator(sat_ice_key_, tag_next_);

  //    and at the current time, where it is a copy evaluator
  S_->Require<CompositeVector,CompositeVectorSpace>(sat_key_, tag_current_, name_);
  // RequireEvaluatorPrimary(sat_key_, tag_current_, *S_);
  S_->Require<CompositeVector,CompositeVectorSpace>(sat_ice_key_, tag_current_, name_);
  // RequireEvaluatorPrimary(sat_ice_key_, tag_current_, *S_);

  // -- the rel perm evaluator, also with the same underlying WRM.
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::BOUNDARY_FACE;
  names2[0] = "cell";
  names2[1] = "boundary_face";
  S_->Require<CompositeVector,CompositeVectorSpace>(coef_key_, tag_next_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponents(names2, locations2, num_dofs2);
  S_->GetEvaluatorList(coef_key_).set<double>("permeability rescaling", perm_scale_);
  S_->RequireEvaluator(coef_key_, tag_next_);

  // -- get the WRM models
  auto wrm_eval = dynamic_cast<Flow::WRMPermafrostEvaluator*>(&wrm);
  AMANZI_ASSERT(wrm_eval != nullptr);
  wrms_ = wrm_eval->get_WRMs();

  // -- molar density used to infer liquid Darcy velocity from flux
  S_->Require<CompositeVector,CompositeVectorSpace>(molar_dens_key_, tag_next_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(molar_dens_key_, tag_next_);

  // -- liquid mass density for the gravity fluxes
  S_->Require<CompositeVector,CompositeVectorSpace>(mass_dens_key_, tag_next_)
    .SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S_->RequireEvaluator(mass_dens_key_, tag_next_); // simply picks up the molar density one.

}

} // namespace
} // namespace
