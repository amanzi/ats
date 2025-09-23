/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Provides a depth-based profile of root density.
#include "rooting_depth_fraction_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
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


double
RootingDepthFractionEvaluator::computeIntegralRootFunc(double z, double alpha, double beta) const
{
  return -0.5 * (std::exp(-alpha * z) + std::exp(-beta * z));
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
  const Epetra_MultiVector& cv = *S.Get<CompositeVector>(cv_key_, tag).ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv =
    *S.Get<CompositeVector>(surf_cv_key_, tag).ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);

  auto& subsurf_mesh = *S.GetMesh(domain_sub_);
  auto& surf_mesh = *S.GetMesh(domain_surf_);

  for (const auto& lc : land_cover_) {
    auto lc_ids = surf_mesh.getSetEntities(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    for (auto sc : lc_ids) {
      double depth = 0.;
      double total = 0.;

      double at_max = 1.0 + computeIntegralRootFunc(lc.second.rooting_depth_max,
                                                    lc.second.rooting_profile_alpha,
                                                    lc.second.rooting_profile_beta);
      const auto& col_cells = subsurf_mesh.columns.getCells(sc);
      int i = 0;
      for (auto c : col_cells) {
        result_v[0][c] = -computeIntegralRootFunc(
                           depth, lc.second.rooting_profile_alpha, lc.second.rooting_profile_beta) /
                         at_max;
        depth += cv[0][c] / surf_cv[0][sc];
        i++;

        if (depth <= lc.second.rooting_depth_max) {
          if (i == (col_cells.size())) {
            // max depth is bigger than the domain, remainder goes in the bottom-most cell
            result_v[0][c] += computeIntegralRootFunc(lc.second.rooting_depth_max,
                                                      lc.second.rooting_profile_alpha,
                                                      lc.second.rooting_profile_beta) /
                              at_max;
            total += result_v[0][c];
          } else {
            result_v[0][c] += computeIntegralRootFunc(depth,
                                                      lc.second.rooting_profile_alpha,
                                                      lc.second.rooting_profile_beta) /
                              at_max;
            total += result_v[0][c];
          }
        } else {
          // max depth stops in this cell, just compute to there
          result_v[0][c] += computeIntegralRootFunc(lc.second.rooting_depth_max,
                                                    lc.second.rooting_profile_alpha,
                                                    lc.second.rooting_profile_beta) /
                            at_max;

          total += result_v[0][c];
          break;
        }
      }
      AMANZI_ASSERT(std::abs(total - 1.0) < 1.e-8);
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
  if (land_cover_.size() == 0) {
    land_cover_ =
      getLandCover(S.GetModelParameters("land cover types"),
                   { "rooting_depth_max", "rooting_profile_alpha", "rooting_profile_beta" });
  }

  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // Create an unowned factory to check my dependencies.
  CompositeVectorSpace dep_fac_one;
  dep_fac_one.SetMesh(S.GetMesh(domain))
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(cv_key_, tag).Update(dep_fac_one);

  CompositeVectorSpace surf_fac_one;
  surf_fac_one.SetMesh(S.GetMesh(Keys::getDomain(surf_cv_key_)))
    ->SetGhosted(true)
    ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
  S.Require<CompositeVector, CompositeVectorSpace>(surf_cv_key_, tag).Update(surf_fac_one);
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
