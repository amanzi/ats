/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Distributes and downregulates potential transpiration to the rooting zone.
#include "errors.hh"
#include "exceptions.hh"
#include "Units.hh"
#include "transpiration_distribution_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

const std::string TranspirationDistributionEvaluator::eval_type =
  "transpiration distribution, rooting depth";

// Constructor from ParameterList
TranspirationDistributionEvaluator::TranspirationDistributionEvaluator(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
TranspirationDistributionEvaluator::Clone() const
{
  return Teuchos::rcp(new TranspirationDistributionEvaluator(*this));
}


// Initialize by setting up dependencies
void
TranspirationDistributionEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  domain_sub_ = Keys::getDomain(my_keys_.front().first);
  domain_surf_ = Keys::readDomainHint(*plist_, domain_sub_, "domain", "surface");
  Tag tag = my_keys_.front().second;

  if (plist_->isSublist("water limiter function")) {
    Errors::Message msg("TranspirationDistributionEvaluator: \"water limiter function\" is not "
                        "supported on device.");
    Exceptions::amanzi_throw(msg);
  }
  limiter_local_ = plist_->get<bool>("water limiter local", true);

  // - pull Keys from plist
  // dependency: pressure
  f_wp_key_ = Keys::readKey(*plist_, domain_sub_, "plant wilting factor", "plant_wilting_factor");
  dependencies_.insert(KeyTag{ f_wp_key_, tag });

  // dependency: rooting_depth_fraction
  f_root_key_ = Keys::readKey(*plist_, domain_sub_, "root fraction", "root_fraction");
  dependencies_.insert(KeyTag{ f_root_key_, tag });

  // dependency: transpiration
  potential_trans_key_ =
    Keys::readKey(*plist_, domain_surf_, "potential transpiration", "potential_transpiration");
  dependencies_.insert(KeyTag{ potential_trans_key_, tag });

  // dependency: cell volume, surface cell volume
  cv_key_ = Keys::readKey(*plist_, domain_sub_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  surf_cv_key_ = Keys::readKey(*plist_, domain_surf_, "surface cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ surf_cv_key_, tag });

  year_duration_ = plist_->get<double>("year duration", 1.0);
  std::string year_duration_units = plist_->get<std::string>("year duration units", "noleap");

  land_cover_ =
    getLandCoverMap(plist_->sublist("model parameters"), { "leaf_on_doy", "leaf_off_doy" });

  // deal with units
  Amanzi::Utils::Units units;
  bool flag;
  year_duration_ = units.ConvertTime(year_duration_, year_duration_units, "s", flag);
}


void
TranspirationDistributionEvaluator::Evaluate_(const State& S,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // on the subsurface
  auto f_wp = S.Get<CompositeVector>(f_wp_key_, tag).viewComponent("cell", false);
  auto f_root = S.Get<CompositeVector>(f_root_key_, tag).viewComponent("cell", false);
  auto cv = S.Get<CompositeVector>(cv_key_, tag).viewComponent("cell", false);

  // on the surface
  auto potential_trans =
    S.Get<CompositeVector>(potential_trans_key_, tag).viewComponent("cell", false);
  auto surf_cv = S.Get<CompositeVector>(surf_cv_key_, tag).viewComponent("cell", false);

  result[0]->putScalar(0.);
  auto result_v = result[0]->viewComponent("cell", false);

  auto surf_mesh = S.GetMesh(domain_surf_);
  const AmanziMesh::MeshCache& subsurf_mesh = S.GetMesh(domain_sub_)->getCache();
  bool limiter_local(limiter_local_);

  for (const auto& region_lc : land_cover_) {
    if (!TranspirationPeriod_(
          S.get_time(), region_lc.second.leaf_on_doy, region_lc.second.leaf_off_doy))
      continue;

    auto lc_ids = surf_mesh->getSetEntities<MemSpace_kind::DEVICE>(
      region_lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    Kokkos::parallel_for(
      "TranspirationDistributionEvaluator::Evaluate", lc_ids.size(), KOKKOS_LAMBDA(const int i) {
        AmanziMesh::Entity_ID sc = lc_ids(i);
        const auto& col_cells = subsurf_mesh.columns.getCells<MemSpace_kind::DEVICE>(sc);

        double column_total = 0.;
        for (auto c : col_cells) {
          column_total += f_wp(c, 0) * f_root(c, 0) * cv(c, 0);
          result_v(c, 0) = f_wp(c, 0) * f_root(c, 0);
        }

        if (column_total > 0.) {
          double coef = potential_trans(sc, 0) * surf_cv(sc, 0) / column_total;
          for (auto c : col_cells) {
            result_v(c, 0) *= coef;
            if (limiter_local) result_v(c, 0) *= f_wp(c, 0);
          }
        }
      });
  }
}


void
TranspirationDistributionEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  result[0]->putScalar(0.); // nonlocal due to rescaling, not computed
}


void
TranspirationDistributionEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(my_keys_.front().first);

  // Create an unowned factory to check my dependencies.
  // -- first those on the subsurface mesh
  CompositeVectorSpace dep_fac;
  dep_fac.SetMesh(S.GetMesh(domain))->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(f_root_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(f_wp_key_, tag).Update(dep_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(cv_key_, tag).Update(dep_fac);

  // -- next those on the surface mesh
  CompositeVectorSpace surf_fac;
  surf_fac.SetMesh(S.GetMesh(Keys::getDomain(surf_cv_key_)))
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(potential_trans_key_, tag).Update(surf_fac);
  S.Require<CompositeVector, CompositeVectorSpace>(surf_cv_key_, tag).Update(surf_fac);
}


bool
TranspirationDistributionEvaluator::TranspirationPeriod_(double time,
                                                         double leaf_on_doy,
                                                         double leaf_off_doy)
{
  if (leaf_on_doy < 0 || leaf_off_doy < 0) {
    return true; // evergreen
  }

  double time_of_year = fmod(time, year_duration_);
  double leaf_on_time = leaf_on_doy * 86400;
  double leaf_off_time = leaf_off_doy * 86400;

  if (leaf_on_time < leaf_off_time) {
    // northern hemisphere
    if ((leaf_on_time <= time_of_year) && (time_of_year < leaf_off_time)) {
      //summer
      return true;
    } else {
      return false;
    }
  } else {
    // southern hemisphere
    if ((leaf_off_time <= time_of_year) && (time_of_year < leaf_on_time)) {
      // southern hemisphere summer
      return true;
    } else {
      return false;
    }
  }
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
