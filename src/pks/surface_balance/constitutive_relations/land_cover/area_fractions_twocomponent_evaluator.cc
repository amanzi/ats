/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "area_fractions_twocomponent_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
AreaFractionsTwoComponentEvaluator::AreaFractionsTwoComponentEvaluator(
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
      "AreaFractionsTwoComponentEvaluator: Minimum fractional area should be > 0.");
    Exceptions::amanzi_throw(message);
  }

  // get domain names
  domain_ = Keys::getDomain(my_keys_.front().first); // surface
  domain_snow_ = Keys::readDomainHint(plist_, domain_, "surface", "snow");
  auto tag = my_keys_.front().second;

  // get dependencies
  snow_depth_key_ = Keys::readKey(plist_, domain_snow_, "snow depth", "depth");
  dependencies_.insert(KeyTag{ snow_depth_key_, tag });
}


void
AreaFractionsTwoComponentEvaluator::Evaluate_(const State& S,
                                              const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  auto mesh = result[0]->Mesh();
  auto& res = *result[0]->ViewComponent("cell", false);
  const auto& sd = *S.Get<CompositeVector>(snow_depth_key_, tag).ViewComponent("cell", false);

  for (const auto& lc : land_cover_) {
    auto lc_ids = mesh->getSetEntities(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (auto c : lc_ids) {
      // calculate area of land
      if (sd[0][c] >= lc.second.snow_transition_depth) {
        res[1][c] = 1.;
      } else if (sd[0][c] <= 0.) {
        res[1][c] = 0.;
      } else {
        res[1][c] = sd[0][c] / lc.second.snow_transition_depth;
      }

      // if any area is less than eps, give to other
      if (res[1][c] < min_area_) {
        res[1][c] = 0.;
      } else if (res[1][c] > (1 - min_area_)) {
        res[1][c] = 1.;
      }
      res[0][c] = 1 - res[1][c];
    }
  }

  // debugging for bad input files
  int nerr = 0;
  for (int c = 0; c != res.MyLength(); ++c) {
    if (std::abs(1 - res[0][c] - res[1][c]) > 1e-10) nerr++;
  }
  int nerr_global = 0;
  mesh->getComm()->SumAll(&nerr, &nerr_global, 1);
  if (nerr_global > 0) {
    Errors::Message msg("AreaFractionsTwoComponent: land cover types do not cover the mesh.");
    Exceptions::amanzi_throw(msg);
  }
}

void
AreaFractionsTwoComponentEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(0.);
  // Errors::Message msg("NotImplemented: AreaFractionsTwoComponentEvaluator currently does not provide derivatives.");
  // Exceptions::amanzi_throw(msg);
}


// custom EC used to set subfield names
void
AreaFractionsTwoComponentEvaluator::EnsureCompatibility_Structure_(State& S)
{
  S.GetRecordSetW(my_keys_.front().first).set_subfieldnames({ "bare_or_water", "snow" });
}


void
AreaFractionsTwoComponentEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (land_cover_.size() == 0)
    land_cover_ = getLandCover(S.ICList().sublist("land cover types"), { "snow_transition_depth" });

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
