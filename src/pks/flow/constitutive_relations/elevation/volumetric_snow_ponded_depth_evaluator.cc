/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

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

#include "MeshAlgorithms.hh"
#include "volumetric_snow_ponded_depth_evaluator.hh"
#include "subgrid_microtopography.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {


VolumetricSnowPondedDepthEvaluator::VolumetricSnowPondedDepthEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  my_keys_.clear(); // clear to push back in order
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
  vol_pd_key_ =
    Keys::readKey(plist, domain_surf_, "volumetric ponded depth", "volumetric_ponded_depth");
  my_keys_.emplace_back(KeyTag{ vol_pd_key_, tag });
  vol_sd_key_ = Keys::readKey(plist, domain_snow_, "volumetric snow depth", "volumetric_depth");
  my_keys_.emplace_back(KeyTag{ vol_sd_key_, tag });

  // dependencies
  pd_key_ = Keys::readKey(plist_, domain_surf_, "ponded depth key", "ponded_depth");
  dependencies_.insert(KeyTag{ pd_key_, tag });
  sd_key_ = Keys::readKey(plist_, domain_snow_, "snow depth key", "depth");
  dependencies_.insert(KeyTag{ sd_key_, tag });

  delta_max_key_ =
    Keys::readKey(plist_, domain_surf_, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(KeyTag{ delta_max_key_, tag });

  delta_ex_key_ = Keys::readKey(plist_, domain_surf_, "excluded volume", "excluded_volume");
  dependencies_.insert(KeyTag{ delta_ex_key_, tag });
}


void
VolumetricSnowPondedDepthEvaluator::Evaluate_(const State& S,
                                              const std::vector<CompositeVector*>& results)
{
  // NOTE, we can only differentiate with respect to quantities that exist on
  // all entities, not just cell entities.
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pd_v = S.GetPtr<CompositeVector>(pd_key_, tag);
  Teuchos::RCP<const CompositeVector> sd_v = S.GetPtr<CompositeVector>(sd_key_, tag);
  Teuchos::RCP<const CompositeVector> del_max_v = S.GetPtr<CompositeVector>(delta_max_key_, tag);
  Teuchos::RCP<const CompositeVector> del_ex_v = S.GetPtr<CompositeVector>(delta_ex_key_, tag);
  const AmanziMesh::Mesh& mesh = *results[0]->Mesh();

  for (const auto& comp : *results[0]) {
    AMANZI_ASSERT(comp == "cell" || comp == "boundary_face");
    bool is_internal_comp = comp == "boundary_face";
    Key internal_comp = is_internal_comp ? "cell" : comp;

    const auto& pd = *pd_v->ViewComponent(comp, false);
    const auto& del_max = *del_max_v->ViewComponent(internal_comp, false);
    const auto& del_ex = *del_ex_v->ViewComponent(internal_comp, false);
    auto& vpd = *results[0]->ViewComponent(comp, false);

    // note, snow depth does not have a boundary face component?
    if (comp == "boundary_face") {
      AMANZI_ASSERT(!results[1]->HasComponent(comp));

      int ncomp = vpd.MyLength();
      for (int i = 0; i != ncomp; ++i) {
        int ii = is_internal_comp ? AmanziMesh::getBoundaryFaceInternalCell(mesh, i) : i;
        AMANZI_ASSERT(Microtopography::validParameters(del_max[0][ii], del_ex[0][ii]));

        double pdc = std::max(0., pd[0][i]);
        vpd[0][i] = Microtopography::volumetricDepth(pdc, del_max[0][ii], del_ex[0][ii]);
      }

    } else if (comp == "cell") {
      auto& vsd = *results[1]->ViewComponent(comp, false);
      const auto& sd = *sd_v->ViewComponent(comp, false);

      int ncomp = vpd.MyLength();
      for (int i = 0; i != ncomp; ++i) {
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
VolumetricSnowPondedDepthEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{
  Tag tag = my_keys_.front().second;
  // NOTE, we can only differentiate with respect to quantities that exist on
  // all entities, not just cell entities.
  //
  // Not differentiating boundary faces... hopefully that is ok.
  const auto& pd = *S.GetPtr<CompositeVector>(pd_key_, tag)->ViewComponent("cell", false);
  const auto& sd = *S.GetPtr<CompositeVector>(sd_key_, tag)->ViewComponent("cell", false);
  const auto& del_max =
    *S.GetPtr<CompositeVector>(delta_max_key_, tag)->ViewComponent("cell", false);
  const auto& del_ex = *S.GetPtr<CompositeVector>(delta_ex_key_, tag)->ViewComponent("cell", false);
  auto& vpd = *results[0]->ViewComponent("cell", false);
  auto& vsd = *results[1]->ViewComponent("cell", false);

  if (wrt_key == pd_key_) {
    for (int c = 0; c != vpd.MyLength(); ++c) {
      double sdc = std::max(0., sd[0][c]);
      double pdc = std::max(0., pd[0][c]);
      vpd[0][c] = Microtopography::dVolumetricDepth_dDepth(pdc, del_max[0][c], del_ex[0][c]);
      vsd[0][c] = Microtopography::dVolumetricDepth_dDepth(pdc + sdc, del_max[0][c], del_ex[0][c]);
    }
  } else if (wrt_key == sd_key_) {
    vpd.PutScalar(0.);
    for (int c = 0; c != vpd.MyLength(); ++c) {
      double sdc = std::max(0., sd[0][c]);
      double pdc = std::max(0., pd[0][c]);
      vsd[0][c] = Microtopography::dVolumetricDepth_dDepth(pdc + sdc, del_max[0][c], del_ex[0][c]);
    }
  } else {
    Errors::Message msg("VolumetricSnowPondedDepthEvaluator: Not Implemented: no derivatives "
                        "implemented other than depths.");
    Exceptions::amanzi_throw(msg);
  }
}

void
VolumetricSnowPondedDepthEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  // Ensure my field exists.  Requirements should be already set.
  const auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                        my_keys_.front().second);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac.Mesh() != Teuchos::null) {
    // Create an unowned factory to check my dependencies.
    Teuchos::RCP<CompositeVectorSpace> dep_fac = Teuchos::rcp(new CompositeVectorSpace(my_fac));
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
    for (const auto& key_tag : dependencies_) {
      if (key_tag.first == delta_max_key_ || key_tag.first == delta_ex_key_ ||
          key_tag.first == sd_key_) {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*no_bf_dep_fac);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
          .Update(*dep_fac);
      }
    }
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
