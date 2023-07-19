/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates albedos and emissivities in a three-region subgrid model.
#include "Key.hh"
#include "albedo_threecomponent_evaluator.hh"
#include "seb_physics_defs.hh"
#include "seb_physics_funcs.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

AlbedoThreeComponentEvaluator::AlbedoThreeComponentEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  // determine the domain
  Key akey = my_keys_.front().first;
  domain_ = Keys::getDomain(akey);
  akey = Keys::getVarName(akey);
  domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");
  auto tag = my_keys_.front().second;

  // my keys
  // -- sources
  my_keys_.clear();
  albedo_key_ = Keys::in(akey, "albedo") ? akey : "albedos";
  albedo_key_ = Keys::readKey(plist, domain_, "albedos", albedo_key_);
  my_keys_.emplace_back(KeyTag{ albedo_key_, tag });

  emissivity_key_ = Keys::in(akey, "emissivit") ? akey : "emissivities";
  emissivity_key_ = Keys::readKey(plist, domain_, "emissivities", emissivity_key_);
  my_keys_.emplace_back(KeyTag{ emissivity_key_, tag });

  // dependencies
  // -- snow properties
  snow_dens_key_ = Keys::readKey(plist, domain_snow_, "snow density", "density");
  dependencies_.insert(KeyTag{ snow_dens_key_, tag });

  // -- skin properties
  unfrozen_fraction_key_ = Keys::readKey(plist, domain_, "unfrozen fraction", "unfrozen_fraction");
  dependencies_.insert(KeyTag{ unfrozen_fraction_key_, tag });

  // parameters
  a_ice_ = plist_.get<double>("albedo ice [-]", 0.44);
  a_water_ = plist_.get<double>("albedo water [-]", 0.1168);

  e_ice_ = plist_.get<double>("emissivity ice [-]", 0.98);
  e_water_ = plist_.get<double>("emissivity water [-]", 0.995);
  e_snow_ = plist_.get<double>("emissivity ground surface [-]", 0.98);
}

// Required methods from EvaluatorSecondaryMonotypeCV
void
AlbedoThreeComponentEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& results)
{
  auto mesh = S.GetMesh(domain_);
  auto tag = my_keys_.front().second;

  // collect dependencies
  const auto& snow_dens = *S.Get<CompositeVector>(snow_dens_key_, tag).ViewComponent("cell", false);
  const auto& unfrozen_fraction =
    *S.Get<CompositeVector>(unfrozen_fraction_key_, tag).ViewComponent("cell", false);

  // collect output vecs
  auto& albedo = *results[0]->ViewComponent("cell", false);
  auto& emissivity = *results[1]->ViewComponent("cell", false);

  emissivity(2)->PutScalar(e_snow_);

  for (const auto& lc : land_cover_) {
    auto lc_ids = mesh->getSetEntities(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (auto c : lc_ids) {
      // albedo of the snow
      albedo[2][c] = Relations::CalcAlbedoSnow(snow_dens[0][c]);

      // a and e of water
      albedo[1][c] = unfrozen_fraction[0][c] * a_water_ + (1 - unfrozen_fraction[0][c]) * a_ice_;
      emissivity[1][c] =
        unfrozen_fraction[0][c] * e_water_ + (1 - unfrozen_fraction[0][c]) * e_ice_;
      // a and e of soil
      albedo[0][c] = lc.second.albedo_ground;
      emissivity[0][c] = lc.second.emissivity_ground;
    }
  }
}

void
AlbedoThreeComponentEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& results)
{}


// custom EC used to set subfield names
void
AlbedoThreeComponentEvaluator::EnsureCompatibility_Structure_(State& S)
{
  S.GetRecordSetW(my_keys_.front().first).set_subfieldnames({ "bare", "water", "snow" });
}


void
AlbedoThreeComponentEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  // new state!
  if (land_cover_.size() == 0)
    land_cover_ = getLandCover(S.ICList().sublist("land cover types"),
                               { "albedo_ground", "emissivity_ground" });

  for (auto dep : dependencies_) {
    auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second);
    if (Keys::getDomain(dep.first) == domain_snow_) {
      fac.SetMesh(S.GetMesh(domain_snow_))->SetGhosted()->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    } else {
      fac.SetMesh(S.GetMesh(domain_))->SetGhosted()->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
