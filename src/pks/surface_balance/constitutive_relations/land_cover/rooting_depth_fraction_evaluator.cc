/*
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
RootingDepthFractionEvaluator::RootingDepthFractionEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
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
  domain_sub_ = Keys::getDomain(my_key_);
  domain_surf_ = Keys::readDomainHint(plist_, domain_sub_, "domain", "surface");

  // - pull Keys from plist
  // dependency: depth
  z_key_ = Keys::readKey(plist_, domain_sub_, "depth", "depth");
  dependencies_.insert(z_key_);

  // cell volume, surface area
  cv_key_ = Keys::readKey(plist_, domain_sub_, "cell volume", "cell_volume");
  dependencies_.insert(cv_key_);

  surf_cv_key_ = Keys::readKey(plist_, domain_surf_, "surface cell volume", "cell_volume");
  dependencies_.insert(surf_cv_key_);
}


void
RootingDepthFractionEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  const Epetra_MultiVector& z = *S->GetFieldData(z_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& cv = *S->GetFieldData(cv_key_)->ViewComponent("cell", false);
  const Epetra_MultiVector& surf_cv = *S->GetFieldData(surf_cv_key_)->ViewComponent("cell", false);
  Epetra_MultiVector& result_v = *result->ViewComponent("cell", false);

  auto& subsurf_mesh = *S->GetMesh(domain_sub_);
  auto& surf_mesh = *S->GetMesh(domain_surf_);

  for (const auto& region_model : models_) {
    AmanziMesh::Entity_ID_List lc_ids;
    surf_mesh.get_set_entities(region_model.first, AmanziMesh::Entity_kind::CELL,
                           AmanziMesh::Parallel_type::OWNED, &lc_ids);

    for (int sc : lc_ids) {
      double column_total = 0.;
      double f_root_total = 0.;
      for (auto c : subsurf_mesh.cells_of_column(sc)) {
        result_v[0][c] = region_model.second->RootingDepthFraction(z[0][c]);
        column_total += result_v[0][c] * cv[0][c];
      }

      // normalize to 1 over the column
      if (column_total > 0) {
        for (auto c : subsurf_mesh.cells_of_column(sc)) {
          result_v[0][c] = result_v[0][c] * 1.0 * surf_cv[0][sc] / column_total;
        }
      }
    }
  }
}


void
RootingDepthFractionEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  // this should only change if the mesh deforms.  don't do that!
  result->PutScalar(0.);
  AMANZI_ASSERT(0);
}


void
RootingDepthFractionEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  if (models_.size() == 0) {
    land_cover_ = getLandCover(S->ICList().sublist("land cover types"),
            {"rooting_profile_alpha", "rooting_profile_beta", "rooting_depth_max"});
    for (const auto& lc : land_cover_) {
      models_[lc.first] = Teuchos::rcp(new RootingDepthFractionModel(lc.second));
    }
  }

  // Ensure my field exists.  Requirements should be already set.
  AMANZI_ASSERT(!my_key_.empty());
  auto my_fac = S->RequireField(my_key_, my_key_);

  // check plist for vis or checkpointing control
  bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
  S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
  bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
  S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

  Key domain = Keys::getDomain(my_key_);
  my_fac->SetMesh(S->GetMesh(domain))
      ->SetComponent("cell", AmanziMesh::CELL, 1)
      ->SetGhosted(true);

  // Create an unowned factory to check my dependencies.
  CompositeVectorSpace dep_fac_one;
  dep_fac_one.SetMesh(my_fac->Mesh())
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(z_key_)->Update(dep_fac_one);
  S->RequireField(cv_key_)->Update(dep_fac_one);

  CompositeVectorSpace surf_fac_one;
  surf_fac_one.SetMesh(S->GetMesh(Keys::getDomain(surf_cv_key_)))
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireField(surf_cv_key_)->Update(surf_fac_one);

  // Recurse into the tree to propagate info to leaves.
  for (KeySet::const_iterator key=dependencies_.begin();
       key!=dependencies_.end(); ++key) {
    S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
  }
}


} //namespace
} //namespace
} //namespace
