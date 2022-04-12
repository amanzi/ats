/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Determine the volumetric ponded depth from ponded depth and subgrid parameters.
/*!

* `"maximum relief key`" ``[string]`` **DOMAIN-maximum_relief**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"ponded depth key`" ``[string]`` **DOMAIN-ponded_depth**
         The true height of the water surface.

*/

#include "Mesh_Algorithms.hh"
#include "volumetric_ponded_depth_evaluator.hh"
#include "subgrid_microtopography.hh"

namespace Amanzi {
namespace Flow {

VolumetricPondedDepthEvaluator::VolumetricPondedDepthEvaluator(Teuchos::ParameterList& plist) :
     SecondaryVariableFieldEvaluator(plist)
{
  Key domain = Keys::getDomain(my_key_);

  // dependencies
  pd_key_ = Keys::readKey(plist_, domain, "ponded depth key", "ponded_depth");
  dependencies_.insert(pd_key_);

  delta_max_key_ = Keys::readKey(plist_, domain, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(delta_max_key_);

  delta_ex_key_ = Keys::readKey(plist_, domain, "excluded volume", "excluded_volume");
  dependencies_.insert(delta_ex_key_);
}


void
VolumetricPondedDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  // NOTE, we can only differentiate with respect to quantities that exist on
  // all entities, not just cell entities.
  Teuchos::RCP<const CompositeVector> pd_v = S->GetFieldData(pd_key_);
  Teuchos::RCP<const CompositeVector> del_max_v = S->GetFieldData(delta_max_key_);
  Teuchos::RCP<const CompositeVector> del_ex_v = S->GetFieldData(delta_ex_key_);
  const AmanziMesh::Mesh& mesh = *result->Mesh();

  for (const auto& comp : *result) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const auto& pd = *pd_v->ViewComponent(comp,false);
    const auto& del_max = *del_max_v->ViewComponent(internal_comp,false);
    const auto& del_ex = *del_ex_v->ViewComponent(internal_comp,false);
    auto& res = *result->ViewComponent(comp,false);

    int ncomp = result->size(comp, false);
    for (int i=0; i!=ncomp; ++i) {
      int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
      AMANZI_ASSERT(Microtopography::validParameters(del_max[0][ii], del_ex[0][ii]));
      res[0][i] = Microtopography::volumetricDepth(pd[0][i], del_max[0][ii], del_ex[0][ii]);
    }
  }
}


void
VolumetricPondedDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> pd_v = S->GetFieldData(pd_key_);
  Teuchos::RCP<const CompositeVector> del_max_v = S->GetFieldData(delta_max_key_);
  Teuchos::RCP<const CompositeVector> del_ex_v = S->GetFieldData(delta_ex_key_);
  const AmanziMesh::Mesh& mesh = *result->Mesh();

  if (wrt_key == pd_key_) {
    for (const auto& comp : *result) {
      AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
      bool is_internal_comp = comp == "boundary_face";
      Key internal_comp = is_internal_comp ? "cell" : comp;

      const auto& pd = *pd_v->ViewComponent(comp,false);
      const auto& del_max = *del_max_v->ViewComponent(internal_comp,false);
      const auto& del_ex = *del_ex_v->ViewComponent(internal_comp,false);
      auto& res = *result->ViewComponent(comp,false);

      int ncomp = result->size(comp, false);
      for (int i=0; i!=ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        res[0][i] = Microtopography::dVolumetricDepth_dDepth(pd[0][i], del_max[0][ii], del_ex[0][ii]);
        res[0][i] = std::max(res[0][i],0.001);
      }
    }

  } else {
    Errors::Message msg("VolumetricPondedDepthEvaluator: Not Implemented: no derivatives implemented other than ponded depth.");
    Exceptions::amanzi_throw(msg);
  }
}

void
VolumetricPondedDepthEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(my_key_ != std::string(""));
  Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac->Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    Teuchos::RCP<CompositeVectorSpace> dep_fac =
        Teuchos::rcp(new CompositeVectorSpace(*my_fac));
    dep_fac->SetOwned(false);

    Teuchos::RCP<CompositeVectorSpace> no_bf_dep_fac;
    if (dep_fac->HasComponent("boundary_face")) {
      no_bf_dep_fac = Teuchos::rcp(new CompositeVectorSpace());
      no_bf_dep_fac->SetMesh(dep_fac->Mesh())
        ->SetGhosted(true)
        ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    } else {
      no_bf_dep_fac = dep_fac;
    }

    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto& key : dependencies_) {
      if (key == my_key_) {
        Errors::Message msg;
        msg << "Evaluator for key \"" << my_key_ << "\" depends upon itself.";
        Exceptions::amanzi_throw(msg);
      }

      if (key == delta_max_key_ || key == delta_ex_key_) {
        Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
        fac->Update(*no_bf_dep_fac);
      } else {
        Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);
        fac->Update(*dep_fac);
      }
    }

    // Recurse into the tree to propagate info to leaves.
    for (const auto& key : dependencies_) {
      S->RequireFieldEvaluator(key)->EnsureCompatibility(S);
    }
  }
}



} //namespace
} //namespace
