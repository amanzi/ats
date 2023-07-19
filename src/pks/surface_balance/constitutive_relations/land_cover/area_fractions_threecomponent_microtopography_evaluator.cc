/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A subgrid model for determining the area fraction of land, water, and snow within a grid cell with subgrid microtopography.
#include "subgrid_microtopography.hh"
#include "area_fractions_threecomponent_microtopography_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
AreaFractionsThreeComponentMicrotopographyEvaluator::
  AreaFractionsThreeComponentMicrotopographyEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  //
  // NOTE: this evaluator simplifies the situation by assuming constant
  // density.  This make it so that ice and water see the same geometry per
  // unit pressure, which isn't quite true thanks to density differences.
  // However, we hypothesize that these differences, on the surface (unlike in
  // the subsurface) really don't matter much. --etc
  snow_subgrid_transition_ = plist_.get<double>("snow transition height [m]", 0.02);
  min_area_ = plist_.get<double>("minimum fractional area [-]", 1.e-5);
  if (min_area_ < 0.) {
    Errors::Message message(
      "AreaFractionsTwoComponentEvaluator: Minimum fractional area should be >= 0.");
    Exceptions::amanzi_throw(message);
  }

  // get domain names
  domain_ = Keys::getDomain(my_keys_.front().first);
  domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");
  auto tag = my_keys_.front().second;

  // get dependencies
  delta_max_key_ =
    Keys::readKey(plist_, domain_, "microtopographic relief", "microtopographic_relief");
  dependencies_.insert(KeyTag{ delta_max_key_, tag });

  delta_ex_key_ = Keys::readKey(plist_, domain_, "excluded volume", "excluded_volume");
  dependencies_.insert(KeyTag{ delta_ex_key_, tag });

  ponded_depth_key_ = Keys::readKey(plist_, domain_, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ ponded_depth_key_, tag });

  snow_depth_key_ = Keys::readKey(plist_, domain_snow_, "snow depth", "depth");
  dependencies_.insert(KeyTag{ snow_depth_key_, tag });

  vol_snow_depth_key_ =
    Keys::readKey(plist_, domain_snow_, "volumetric snow depth", "volumetric_depth");
  dependencies_.insert(KeyTag{ vol_snow_depth_key_, tag });
}


void
AreaFractionsThreeComponentMicrotopographyEvaluator::Evaluate_(
  const State& S,
  const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  Epetra_MultiVector& res = *result[0]->ViewComponent("cell", false);

  const Epetra_MultiVector& pd =
    *S.Get<CompositeVector>(ponded_depth_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& sd =
    *S.Get<CompositeVector>(snow_depth_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& vsd =
    *S.Get<CompositeVector>(vol_snow_depth_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& del_max =
    *S.Get<CompositeVector>(delta_max_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& del_ex =
    *S.Get<CompositeVector>(delta_ex_key_, tag).ViewComponent("cell", false);

  for (int c = 0; c != res.MyLength(); ++c) {
    // calculate area of land
    AMANZI_ASSERT(Flow::Microtopography::validParameters(del_max[0][c], del_ex[0][c]));
    double liquid_water_area =
      Flow::Microtopography::dVolumetricDepth_dDepth(pd[0][c], del_max[0][c], del_ex[0][c]);
    double wet_area = Flow::Microtopography::dVolumetricDepth_dDepth(
      pd[0][c] + std::max(sd[0][c], 0.0), del_max[0][c], del_ex[0][c]);

    // now partition the wet area into snow and water
    if (vsd[0][c] >= wet_area * snow_subgrid_transition_) {
      res[2][c] = wet_area;
      res[1][c] = 0.;
      res[0][c] = 1 - wet_area;
    } else {
      res[2][c] = vsd[0][c] / snow_subgrid_transition_;

      // how much of the remainder goes to water?
      res[1][c] = std::min(wet_area - res[2][c], liquid_water_area);
      res[0][c] = 1 - res[1][c] - res[2][c];
    }

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

    AMANZI_ASSERT(std::abs(res[0][c] + res[1][c] + res[2][c] - 1.0) < 1.e-6);
    AMANZI_ASSERT(-1.e-10 <= res[0][c] && res[0][c] <= 1. + 1.e-10);
    AMANZI_ASSERT(-1.e-10 <= res[1][c] && res[1][c] <= 1. + 1.e-10);
    AMANZI_ASSERT(-1.e-10 <= res[2][c] && res[2][c] <= 1. + 1.e-10);

    res[0][c] = std::min(std::max(0., res[0][c]), 1.);
    res[1][c] = std::min(std::max(0., res[1][c]), 1.);
    res[2][c] = std::min(std::max(0., res[2][c]), 1.);
  }
}


// custom EC used to set subfield names
void
AreaFractionsThreeComponentMicrotopographyEvaluator::EnsureCompatibility_Structure_(State& S)
{
  S.GetRecordSetW(my_keys_.front().first).set_subfieldnames({ "bare", "water", "snow" });
}


void
AreaFractionsThreeComponentMicrotopographyEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  auto tag = my_keys_.front().second;
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
