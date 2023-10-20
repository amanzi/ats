/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Provides a depth-based profile of root density.
#include "rooting_depth_fraction_evaluator.hh"
#include "rooting_depth_fraction_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
RootingDepthFractionEvaluator::RootingDepthFractionEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
RootingDepthFractionEvaluator::Clone() const
{
  return Teuchos::rcp(new RootingDepthFractionEvaluator(*this));
}


// Initialize by setting up dependencies
void
RootingDepthFractionEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  domain_sub_ = Keys::getDomain(my_keys_.front().first);
  domain_surf_ = Keys::readDomainHint(plist_, domain_sub_, "domain", "surface");
  Tag tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: depth
  z_key_ = Keys::readKey(plist_, domain_sub_, "depth", "depth");
  dependencies_.insert(KeyTag{ z_key_, tag });

  // cell volume, surface area
  cv_key_ = Keys::readKey(plist_, domain_sub_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  surf_cv_key_ = Keys::readKey(plist_, domain_surf_, "surface cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ surf_cv_key_, tag });
}


void
RootingDepthFractionEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& z = *S.Get<CompositeVector>(z_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv =
    *S.Get<CompositeVector>(surf_cv_key_, tag).ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  auto& subsurf_mesh = *S.GetMesh(domain_sub_);
  auto& surf_mesh = *S.GetMesh(domain_surf_);

  for (const auto& region_model : models_) {
    auto lc_ids = surf_mesh.getSetEntities(
      region_model.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (int sc : lc_ids) {
      double column_total = 0.;
      double f_root_total = 0.;
      for (auto c : subsurf_mesh.columns.getCells(sc)) {
        result_v[0][c] = region_model.second->RootingDepthFraction(z[0][c]);
        column_total += result_v[0][c] * cv[0][c];
      }

      // normalize to 1 over the column
      if (column_total > 0) {
        for (auto c : subsurf_mesh.columns.getCells(sc)) {
          result_v[0][c] = result_v[0][c] * 1.0 * surf_cv[0][sc] / column_total;
        }
      }
    }
  }
}


void
RootingDepthFractionEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  // this should only change if the mesh deforms.  don't do that!
  result[0]->PutScalar(0.);
}


void
RootingDepthFractionEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (models_.size() == 0) {
    land_cover_ =
      getLandCover(S.ICList().sublist("land cover types"),
                   { "rooting_profile_alpha", "rooting_profile_beta", "rooting_depth_max" });
    for (const auto& lc : land_cover_) {
      models_[lc.first] = Teuchos::rcp(new RootingDepthFractionModel(lc.second));
    }
  }

  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // Create an unowned factory to check my dependencies.
  CompositeVectorSpace dep_fac_one;
  dep_fac_one.SetMesh(S.GetMesh(domain))
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);

  S.Require<CompositeVector, CompositeVectorSpace>(z_key_, tag).Update(dep_fac_one);
  S.Require<CompositeVector, CompositeVectorSpace>(cv_key_, tag).Update(dep_fac_one);

  CompositeVectorSpace surf_fac_one;
  surf_fac_one.SetMesh(S.GetMesh(Keys::getDomain(surf_cv_key_)))
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(surf_cv_key_, tag).Update(surf_fac_one);
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace Amanzi
