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
namespace SurfaceBalance {
namespace Relations {

const std::string RootingDepthFractionEvaluator::eval_type = "root fraction";

// Constructor from ParameterList
RootingDepthFractionEvaluator::RootingDepthFractionEvaluator(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
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
  domain_surf_ = Keys::readDomainHint(*plist_, domain_sub_, "domain", "surface");
  Tag tag = my_keys_.front().second;

  cv_key_ = Keys::readKey(*plist_, domain_sub_, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  surf_cv_key_ = Keys::readKey(*plist_, domain_surf_, "surface cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ surf_cv_key_, tag });

  land_cover_ =
    getLandCoverMap(plist_->sublist("model parameters"),
                    { "rooting_depth_max", "rooting_profile_alpha", "rooting_profile_beta" });
}


void
RootingDepthFractionEvaluator::Evaluate_(const State& S,
                                         const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  auto cv = S.Get<CompositeVector>(cv_key_, tag).viewComponent("cell", false);
  auto surf_cv = S.Get<CompositeVector>(surf_cv_key_, tag).viewComponent("cell", false);
  auto result_v = result[0]->viewComponent("cell", false);

  auto surf_mesh_host = S.GetMesh(domain_surf_);
  const AmanziMesh::MeshCache& subsurf_mesh = S.GetMesh(domain_sub_)->getCache();

  for (const auto& lc : land_cover_) {
    const LandCover& lc_pars = lc.second;
    auto lc_ids = surf_mesh_host->getSetEntities<MemSpace_kind::DEVICE>(
      lc.first, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    Kokkos::parallel_for(
      "RootingDepthFractionEvaluator::Evaluate", lc_ids.extent(0), KOKKOS_LAMBDA(const int j) {
        AmanziMesh::Entity_ID sc = lc_ids(j);
        double depth = 0.;
        double total = 0.;

        double at_max = 1.0 + Impl::computeIntegralRootFunc(lc_pars.rooting_depth_max,
                                                      lc_pars.rooting_profile_alpha,
                                                      lc_pars.rooting_profile_beta);
        const auto& col_cells = subsurf_mesh.columns.getCells<MemSpace_kind::DEVICE>(sc);
        int i = 0;
        for (auto c : col_cells) {
          result_v(c, 0) = -Impl::computeIntegralRootFunc(
                             depth, lc_pars.rooting_profile_alpha, lc_pars.rooting_profile_beta) /
                           at_max;
          depth += cv(c, 0) / surf_cv(sc, 0);
          i++;

          if (depth <= lc_pars.rooting_depth_max) {
            if (i == (col_cells.size())) {
              // max depth is bigger than the domain, remainder goes in the bottom-most cell
              result_v(c, 0) += Impl::computeIntegralRootFunc(lc_pars.rooting_depth_max,
                                                        lc_pars.rooting_profile_alpha,
                                                        lc_pars.rooting_profile_beta) /
                                at_max;
              total += result_v(c, 0);
            } else {
              result_v(c, 0) += Impl::computeIntegralRootFunc(depth,
                                                        lc_pars.rooting_profile_alpha,
                                                        lc_pars.rooting_profile_beta) /
                                at_max;
              total += result_v(c, 0);
            }
          } else {
            // max depth stops in this cell, just compute to there
            result_v(c, 0) += Impl::computeIntegralRootFunc(lc_pars.rooting_depth_max,
                                                      lc_pars.rooting_profile_alpha,
                                                      lc_pars.rooting_profile_beta) /
                              at_max;

            total += result_v(c, 0);
            break;
          }
        }
        assert(std::abs(total - 1.0) < 1.e-8);
      });
  }
}


void
RootingDepthFractionEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(false);
}


void
RootingDepthFractionEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
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
} // namespace Amanzi
