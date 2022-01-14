/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
This is the flow component of the Amanzi code. 
License: BSD
Authors: Ethan Coon (ecoon@lanl.gov)
------------------------------------------------------------------------- */


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

#include "permafrost.hh"

namespace Amanzi {
namespace Flow {

// -------------------------------------------------------------
// Create the physical evaluators for water content, water
// retention, rel perm, etc, that are specific to Richards.
// -------------------------------------------------------------
void Permafrost::SetupPhysicalEvaluators_(const Teuchos::Ptr<State>& S) {
  // -- Absolute permeability.
  //       For now, we assume scalar permeability.  This will change.
  S->Require<CompositeVector,CompositeVectorSpace>(perm_key_, Tags::NEXT).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(perm_key_);
  
  S->Require<CompositeVector,CompositeVectorSpace>(conserved_key_, Tags::NEXT).SetMesh(mesh_)->SetGhosted()
    ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(conserved_key_);

  // -- Water retention evaluators, for saturation and rel perm.
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::BOUNDARY_FACE;
  names2[0] = "cell";
  names2[1] = "boundary_face";

  // -- rel perm on cells + boundary faces

  S->Require<CompositeVector,CompositeVectorSpace>(coef_key_, Tags::NEXT).SetMesh(mesh_)->SetGhosted()
      ->AddComponents(names2,locations2,num_dofs2);
 
  // -- This setup is a little funky -- we use four evaluators to capture the physics.
  Teuchos::ParameterList wrm_plist = plist_->sublist("water retention evaluator");
  wrm_plist.set("evaluator name", sat_key_);
  Teuchos::RCP<Flow::WRMPermafrostEvaluator> wrm =
      Teuchos::rcp(new Flow::WRMPermafrostEvaluator(wrm_plist));
  

  if (!S->HasEvaluator(sat_key_)) {
    S->SetEvaluator(sat_key_, wrm);
    S->SetEvaluator(sat_gas_key_, wrm);
    S->SetEvaluator(sat_ice_key_, wrm);
  }

  // -- the rel perm evaluator, also with the same underlying WRM.
  wrm_plist.set("permeability rescaling", perm_scale_);
  wrm_plist.setName(coef_key_);
  wrm_plist.set("evaluator name", coef_key_);
  Teuchos::RCP<Flow::RelPermEvaluator> rel_perm_evaluator =
      Teuchos::rcp(new Flow::RelPermEvaluator(wrm_plist, wrm->get_WRMs()));
  wrms_ = wrm->get_WRMs();
  


  S->SetEvaluator(coef_key_, rel_perm_evaluator);
  
  // -- Liquid density and viscosity for the transmissivity.

  S->Require<CompositeVector,CompositeVectorSpace>(molar_dens_key_, Tags::NEXT).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(molar_dens_key_);

  /* S->Require<CompositeVector,CompositeVectorSpace>("viscosity_liquid", Tags::NEXT).SetMesh(S->GetMesh())->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator("viscosity_liquid");
  */
  // -- liquid mass density for the gravity fluxes
  S->Require<CompositeVector,CompositeVectorSpace>(mass_dens_key_, Tags::NEXT).SetMesh(mesh_)->SetGhosted()
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireEvaluator(mass_dens_key_); // simply picks up the molar density one.

}

} // namespace
} // namespace
