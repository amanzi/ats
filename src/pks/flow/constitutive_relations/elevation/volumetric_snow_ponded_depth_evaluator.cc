/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  License: see $ATS_DIR/COPYRIGHT
  Authors: Ahmad Jan (jana@ornl.gov)
           Ethan Coon (ecoon@ornl.gov)
*/

//! Determine the volumetric ponded and snow depths from ponded depth and snow depth.
/*!

* `"maximum relief key`" ``[string]`` **DOMAIN-maximum_relief**
         The name of del_max, the max microtopography value.
* `"excluded volume key`" ``[string]`` **DOMAIN-excluded_volume**
         The name of del_excluded, the integral of the microtopography.
* `"ponded depth key`" ``[string]`` **DOMAIN-ponded_depth**
         The true height of the water surface.
* `"snow depth key`" ``[string]`` **SNOW_DOMAIN-depth**
         The true height of the snow surface.

*/

#include "Mesh_Algorithms.hh"
#include "volumetric_snow_ponded_depth_evaluator.hh"
#include "subgrid_microtopography.hh"

namespace Amanzi {
namespace Flow {


VolumetricSnowPondedDepthEvaluator::VolumetricSnowPondedDepthEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist)
{
  Key domain = Keys::getDomain(Keys::cleanPListName(plist_.name()));
  Key dtype = Keys::guessDomainType(domain);
  if (dtype == "surface") {
    domain_surf_ = domain;
    domain_snow_ = Keys::readDomainHint(plist_, domain_surf_, "surface", "snow");
  } else if (dtype == "snow") {
    domain_snow_ = domain;
    domain_surf_ = Keys::readDomainHint(plist_, domain_snow_, "snow", "surface");
  } else {
    Errors::Message msg("VolumetricSnowPondedDepthEvaluator: not sure how to interpret domain.");
    Exceptions::amanzi_throw(msg);
  }

  // my keys
  vol_pd_key_ = Keys::readKey(plist, domain_surf_, "volumetric ponded depth", "volumetric_ponded_depth");
  my_keys_.push_back(vol_pd_key_);
  vol_sd_key_ = Keys::readKey(plist, domain_snow_, "volumetric snow depth", "volumetric_depth");
  my_keys_.push_back(vol_sd_key_);

  // dependencies
  pd_key_ = Keys::readKey(plist_, domain_surf_, "ponded depth key", "ponded_depth");
  dependencies_.insert(pd_key_);
  sd_key_ = Keys::readKey(plist_, domain_snow_, "snow depth key", "depth");
  dependencies_.insert(sd_key_);

  delta_max_key_ = Keys::readKey(plist_, domain_surf_, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(delta_max_key_);

  delta_ex_key_ = Keys::readKey(plist_, domain_surf_, "excluded volume", "excluded_volume");
  dependencies_.insert(delta_ex_key_);
}


void
VolumetricSnowPondedDepthEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                             const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  // NOTE, we can only differentiate with respect to quantities that exist on
  // all entities, not just cell entities.
  Teuchos::RCP<const CompositeVector> pd_v = S->GetFieldData(pd_key_);
  Teuchos::RCP<const CompositeVector> sd_v = S->GetFieldData(sd_key_);
  Teuchos::RCP<const CompositeVector> del_max_v = S->GetFieldData(delta_max_key_);
  Teuchos::RCP<const CompositeVector> del_ex_v = S->GetFieldData(delta_ex_key_);
  const AmanziMesh::Mesh& mesh = *results[0]->Mesh();

  for (const auto& comp : *results[0]) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const auto& pd = *pd_v->ViewComponent(comp,false);
    const auto& del_max = *del_max_v->ViewComponent(internal_comp,false);
    const auto& del_ex = *del_ex_v->ViewComponent(internal_comp,false);
    auto& vpd = *results[0]->ViewComponent(comp,false);

    // note, snow depth does not have a boundary face component?
    if (comp == "boundary_face") {
      AMANZI_ASSERT(!results[1]->HasComponent(comp));

      int ncomp = vpd.MyLength();
      for (int i=0; i!=ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        AMANZI_ASSERT(Microtopography::validParameters(del_max[0][ii], del_ex[0][ii]));

        double pdc = std::max(0., pd[0][i]);
        vpd[0][i] = Microtopography::volumetricDepth(pdc, del_max[0][ii], del_ex[0][ii]);
      }

    } else if (comp == "cell") {
      auto& vsd = *results[1]->ViewComponent(comp,false);
      const auto& sd = *sd_v->ViewComponent(comp,false);

      int ncomp = vpd.MyLength();
      for (int i=0; i!=ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        AMANZI_ASSERT(Microtopography::validParameters(del_max[0][ii], del_ex[0][ii]));

        double pdc = std::max(0., pd[0][i]);
        double sdc = std::max(0., sd[0][i]);
        vpd[0][i] = Microtopography::volumetricDepth(pdc, del_max[0][ii], del_ex[0][ii]);
        double vol_tot = Microtopography::volumetricDepth(pdc + sdc, del_max[0][ii], del_ex[0][ii]);
        vsd[0][i] = vol_tot - vpd[0][i];
      }
    }
  }
}


void
VolumetricSnowPondedDepthEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> >& results)
{
  // NOTE, we can only differentiate with respect to quantities that exist on
  // all entities, not just cell entities.
  //
  // Not differentiating boundary faces... hopefully that is ok.
  const auto& pd = *S->GetFieldData(pd_key_)->ViewComponent("cell", false);
  const auto& sd = *S->GetFieldData(sd_key_)->ViewComponent("cell", false);
  const auto& del_max = *S->GetFieldData(delta_max_key_)->ViewComponent("cell", false);
  const auto& del_ex = *S->GetFieldData(delta_ex_key_)->ViewComponent("cell", false);
  auto& vpd = *results[0]->ViewComponent("cell",false);
  auto& vsd = *results[1]->ViewComponent("cell",false);

  if (wrt_key == pd_key_) {
    for (int c=0; c!=vpd.MyLength(); ++c) {
      double sdc = std::max(0., sd[0][c]);
      double pdc = std::max(0., pd[0][c]);
      vpd[0][c] = Microtopography::dVolumetricDepth_dDepth(pdc, del_max[0][c], del_ex[0][c]);
      vsd[0][c] = Microtopography::dVolumetricDepth_dDepth(pdc + sdc, del_max[0][c], del_ex[0][c]);
    }
  } else if (wrt_key == sd_key_) {
    vpd.PutScalar(0.);
    for (int c=0; c!=vpd.MyLength(); ++c) {
      double sdc = std::max(0., sd[0][c]);
      double pdc = std::max(0., pd[0][c]);
      vsd[0][c] = Microtopography::dVolumetricDepth_dDepth(pdc + sdc, del_max[0][c], del_ex[0][c]);
    }
  } else {
    Errors::Message msg("VolumetricSnowPondedDepthEvaluator: Not Implemented: no derivatives implemented other than depths.");
    Exceptions::amanzi_throw(msg);
  }
}

void
VolumetricSnowPondedDepthEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  // require my keys
  auto my_fac = S->RequireField(vol_pd_key_, vol_pd_key_);
  my_fac->SetMesh(S->GetMesh(domain_surf_))
    ->SetGhosted();
  //   ->SetComponent("cell", AmanziMesh::CELL, 1);

  auto my_fac_snow = S->RequireField(vol_sd_key_, vol_sd_key_);
  // my_fac_snow->SetOwned(false); // why was this here?
  my_fac_snow->SetMesh(S->GetMesh(domain_snow_))
    ->SetGhosted();
  //   ->SetComponent("cell", AmanziMesh::CELL, 1);

  // Check plist for vis or checkpointing control.
  bool io_my_key = plist_.get<bool>("visualize", true);
  S->GetField(vol_pd_key_, vol_pd_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
  S->GetField(vol_pd_key_, vol_pd_key_)->set_io_checkpoint(checkpoint_my_key);

  for (auto dep_key : dependencies_) {
    auto fac = S->RequireField(dep_key);
    if (Keys::getDomain(dep_key) == domain_snow_) {
      fac->SetMesh(S->GetMesh(domain_snow_))
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    } else {
      fac->SetMesh(S->GetMesh(domain_surf_))
        ->SetGhosted()
        ->AddComponent("cell", AmanziMesh::CELL, 1);
    }

    // Recurse into the tree to propagate info to leaves.
    S->RequireFieldEvaluator(dep_key)->EnsureCompatibility(S);
  }
}

} //namespace
} //namespace
