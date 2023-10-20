/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Distributes and downregulates potential transpiration to the rooting zone.
#include "Function.hh"
#include "FunctionFactory.hh"
#include "transpiration_distribution_evaluator.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
TranspirationDistributionEvaluator::TranspirationDistributionEvaluator(
  Teuchos::ParameterList& plist)
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
  domain_surf_ = Keys::readDomainHint(plist_, domain_sub_, "domain", "surface");
  Tag tag = my_keys_.front().second;

  limiter_local_ = false;
  if (plist_.isSublist("water limiter function")) {
    Amanzi::FunctionFactory fac;
    limiter_ = Teuchos::rcp(fac.Create(plist_.sublist("water limiter function")));
  } else {
    limiter_local_ = plist_.get<bool>("water limiter local", true);
  }

  // - pull Keys from plist
  // dependency: pressure
  f_wp_key_ = Keys::readKey(plist_, domain_sub_, "plant wilting factor", "plant_wilting_factor");
  dependencies_.insert(KeyTag{ f_wp_key_, tag });

  // dependency: rooting_depth_fraction
  f_root_key_ =
    Keys::readKey(plist_, domain_sub_, "rooting depth fraction", "rooting_depth_fraction");
  dependencies_.insert(KeyTag{ f_root_key_, tag });

  // dependency: transpiration
  potential_trans_key_ =
    Keys::readKey(plist_, domain_surf_, "potential transpiration", "potential_transpiration");
  dependencies_.insert(KeyTag{ potential_trans_key_, tag });

  // dependency: cell volume, surface cell volume
  cv_key_ = Keys::readKey(plist_, domain_sub_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
  surf_cv_key_ = Keys::readKey(plist_, domain_surf_, "surface cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ surf_cv_key_, tag });

  year_duration_ = plist_.get<double>("year duration", 1.0);
  std::string year_duration_units = plist_.get<std::string>("year duration units", "noleap");

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
  const Epetra_MultiVector& f_wp =
    *S.Get<CompositeVector>(f_wp_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& f_root =
    *S.Get<CompositeVector>(f_root_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);

  // on the surface
  const Epetra_MultiVector& potential_trans =
    *S.Get<CompositeVector>(potential_trans_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv =
    *S.Get<CompositeVector>(surf_cv_key_, tag).ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  double p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);

  auto& subsurf_mesh = *S.GetMesh(domain_sub_);
  auto& surf_mesh = *S.GetMesh(domain_surf_);

  result_v.PutScalar(0.);
  for (const auto& region_lc : land_cover_) {
    auto lc_ids = surf_mesh.getSetEntities(
      region_lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    if (TranspirationPeriod_(
          S.get_time(), region_lc.second.leaf_on_doy, region_lc.second.leaf_off_doy)) {
      for (int sc : lc_ids) {
        double column_total = 0.;
        double f_root_total = 0.;
        double f_wp_total = 0.;
        double var_dz = 0.;
        for (auto c : subsurf_mesh.columns.getCells(sc)) {
          column_total += f_wp[0][c] * f_root[0][c] * cv[0][c];
          result_v[0][c] = f_wp[0][c] * f_root[0][c];
          if (f_wp[0][c] * f_root[0][c] > 0) var_dz += cv[0][c];
        }

        if (column_total > 0.) {
          double coef = potential_trans[0][sc] * surf_cv[0][sc] / column_total;
          if (limiter_.get()) {
            auto column_total_vector = std::vector<double>(1, column_total / surf_cv[0][sc]);
            double limiting_factor = (*limiter_)(column_total_vector);
            AMANZI_ASSERT(limiting_factor >= 0.);
            AMANZI_ASSERT(limiting_factor <= 1.);
            coef *= limiting_factor;
          }

          for (auto c : subsurf_mesh.columns.getCells(sc)) {
            result_v[0][c] *= coef;
            if (limiter_local_) { result_v[0][c] *= f_wp[0][c]; }
          }
        }
      }
    }
  }
}


void
TranspirationDistributionEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(
    0.); // this would be a nontrivial calculation, as it is technically nonlocal due to rescaling issues?
}


void
TranspirationDistributionEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Tag tag = my_keys_.front().second;

  // new state!
  if (land_cover_.size() == 0)
    land_cover_ =
      getLandCover(S.ICList().sublist("land cover types"), { "leaf_on_doy", "leaf_off_doy" });

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
