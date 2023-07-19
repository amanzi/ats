/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! A subgrid model for determining the area fraction of land, water, and snow within a grid cell with subgrid microtopography.
#include "area_fractions_threecomponent_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
AreaFractionsThreeComponentEvaluator::AreaFractionsThreeComponentEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  //
  // NOTE: this evaluator simplifies the situation by assuming constant
  // density.  This make it so that ice and water see the same geometry per
  // unit pressure, which isn't quite true thanks to density differences.
  // However, we hypothesize that these differences, on the surface (unlike in
  // the subsurface) really don't matter much. --etc
  min_area_ = plist_.get<double>("minimum fractional area [-]", 1.e-5);
  if (min_area_ <= 0.) {
    Errors::Message message(
      "AreaFractionsThreeComponentEvaluator: Minimum fractional area should be > 0.");
    Exceptions::amanzi_throw(message);
  }

  // get domain names
  domain_ = Keys::getDomain(my_keys_.front().first); // surface
  domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");
  auto tag = my_keys_.front().second;

  snow_depth_key_ = Keys::readKey(plist_, domain_snow_, "snow depth", "depth");
  dependencies_.insert(KeyTag{ snow_depth_key_, tag });
  ponded_depth_key_ = Keys::readKey(plist_, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ ponded_depth_key_, tag });
}


void
AreaFractionsThreeComponentEvaluator::Evaluate_(const State& S,
                                                const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  auto mesh = result[0]->Mesh();
  auto& res = *result[0]->ViewComponent("cell", false);
  const auto& sd = *S.Get<CompositeVector>(snow_depth_key_, tag).ViewComponent("cell", false);
  const auto& pd = *S.Get<CompositeVector>(ponded_depth_key_, tag).ViewComponent("cell", false);

  for (const auto& lc : land_cover_) {
    auto lc_ids = mesh->getSetEntities(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (auto c : lc_ids) {
      // calculate area of land
      if (sd[0][c] >= lc.second.snow_transition_depth) {
        res[1][c] = 0.;
        res[2][c] = 1.;
      } else {
        if (sd[0][c] <= 0.) {
          res[2][c] = 0.;
        } else {
          res[2][c] = sd[0][c] / lc.second.snow_transition_depth;
        }

        // snow preferentially covers water, as both go to low lying areas
        if (pd[0][c] >= lc.second.water_transition_depth) {
          res[1][c] = 1 - res[2][c];
        } else if (pd[0][c] <= 0.) {
          res[1][c] = 0.;
        } else {
          double water_covered = pd[0][c] / lc.second.water_transition_depth;
          if (res[2][c] > water_covered) {
            res[1][c] = 0;
          } else {
            res[1][c] = water_covered - res[2][c];
          }
        }
      }
      res[0][c] = 1 - res[1][c] - res[2][c];

      // if any area is less than eps, give to others
      // if any area fraction is less than eps, give it to the others
      if (res[0][c] > 0 && res[0][c] < min_area_) {
        if (res[1][c] < min_area_) {
          res[2][c] = 1.;
          res[1][c] = 0.;
          res[0][c] = 0.;
        } else {
          res[1][c] += res[0][c] * res[1][c] / (res[1][c] + res[2][c]);
          res[2][c] += res[0][c] * res[2][c] / (res[1][c] + res[2][c]);
          res[0][c] = 0.;
        }
      } else if (res[1][c] > 0 && res[1][c] < min_area_) {
        if (res[2][c] < min_area_) {
          res[0][c] = 1.;
          res[1][c] = 0.;
          res[2][c] = 0.;
        } else {
          res[0][c] += res[1][c] * res[0][c] / (res[0][c] + res[2][c]);
          res[2][c] += res[1][c] * res[2][c] / (res[0][c] + res[2][c]);
          res[1][c] = 0.;
        }
      } else if (res[2][c] > 0 && res[2][c] < min_area_) {
        res[0][c] += res[2][c] * res[0][c] / (res[0][c] + res[1][c]);
        res[1][c] += res[2][c] * res[1][c] / (res[0][c] + res[1][c]);
        res[2][c] = 0.;
      }

      AMANZI_ASSERT(std::abs(res[0][c] + res[1][c] + res[2][c] - 1.0) < 1.e-10);
      AMANZI_ASSERT(-1.e-10 <= res[0][c] && res[0][c] <= 1. + 1.e-10);
      AMANZI_ASSERT(-1.e-10 <= res[1][c] && res[1][c] <= 1. + 1.e-10);
      AMANZI_ASSERT(-1.e-10 <= res[2][c] && res[1][c] <= 1. + 1.e-10);

      res[0][c] = std::min(std::max(0., res[0][c]), 1.);
      res[1][c] = std::min(std::max(0., res[1][c]), 1.);
      res[2][c] = std::min(std::max(0., res[2][c]), 1.);
    }
  }

  // debugging for bad input files
  int nerr = 0;
  for (int c = 0; c != res.MyLength(); ++c) {
    if (std::abs(1 - res[0][c] - res[1][c] - res[2][c]) > 1e-10) nerr++;
  }
  int nerr_global = 0;
  mesh->getComm()->SumAll(&nerr, &nerr_global, 1);
  if (nerr_global > 0) {
    Errors::Message msg("AreaFractionsTwoComponent: land cover types do not cover the mesh.");
    Exceptions::amanzi_throw(msg);
  }
}

// custom EC used to set subfield names
void
AreaFractionsThreeComponentEvaluator::EnsureCompatibility_Structure_(State& S)
{
  S.GetRecordSetW(my_keys_.front().first).set_subfieldnames({ "bare", "water", "snow" });
}


void
AreaFractionsThreeComponentEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (land_cover_.size() == 0)
    land_cover_ = getLandCover(S.ICList().sublist("land cover types"),
                               { "snow_transition_depth", "water_transition_depth" });

  auto tag = my_keys_.front().second;
  for (auto& dep : dependencies_) {
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
