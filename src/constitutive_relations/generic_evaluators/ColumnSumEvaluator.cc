/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ahmad Jan (jana@ornl.gov)
*/
//! Sums a subsurface field vertically only a surface field.

#include "ColumnSumEvaluator.hh"

namespace Amanzi {
namespace Relations {

ColumnSumEvaluator::ColumnSumEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotypeCV(plist)
{
  surf_domain_ = Keys::getDomain(my_keys_.front().first);
  domain_ = Keys::readDomainHint(plist_, surf_domain_, "surface", "subsurface");

  dep_key_ = Keys::readKey(plist_, domain_, "summed",
                           Keys::getKey(domain_, Keys::getVarName(my_keys_.front().first)));
  Tag tag = my_keys_.front().second;
  dependencies_.insert(KeyTag{dep_key_, tag});

  Key pname = dep_key_ + " coefficient";
  coef_ = plist_.get<double>(pname, 1.0);

  // dependency: cell volume, surface cell volume
  if (plist_.get<bool>("include volume factor", true)) {
    cv_key_ = Keys::readKey(plist_, domain_, "cell volume", "cell_volume");
    dependencies_.insert(KeyTag{cv_key_, tag});

    surf_cv_key_ = Keys::readKey(plist_, surf_domain_, "surface cell volume", "cell_volume");
    dependencies_.insert(KeyTag{surf_cv_key_, tag});
  }

  if (plist_.get<bool>("divide by density", true)) {
    molar_dens_key_ = Keys::readKey(plist_, domain_, "molar density", "molar_density_liquid");
    dependencies_.insert(KeyTag{molar_dens_key_, tag});
  }
}


Teuchos::RCP<Evaluator>
ColumnSumEvaluator::Clone() const
{
  return Teuchos::rcp(new ColumnSumEvaluator(*this));
}


void
ColumnSumEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& dep_c = *S.Get<CompositeVector>(dep_key_, tag)
    .ViewComponent("cell", false);

  if (cv_key_ != "") {
    const Epetra_MultiVector& cv = *S.Get<CompositeVector>(cv_key_, tag)
      .ViewComponent("cell", false);
    const Epetra_MultiVector& surf_cv = *S.Get<CompositeVector>(surf_cv_key_, tag)
      .ViewComponent("cell", false);

    if (molar_dens_key_ != "") {
      const Epetra_MultiVector& dens = *S.Get<CompositeVector>(molar_dens_key_, tag)
        .ViewComponent("cell",false);

      auto subsurf_mesh = S.GetMesh(domain_);
      for (int c=0; c!=res_c.MyLength(); ++c) {
        double sum = 0;
        for (auto i : subsurf_mesh->cells_of_column(c)) {
          sum += dep_c[0][i] * cv[0][i] / dens[0][i];
        }
        res_c[0][c] = coef_ * sum / surf_cv[0][c];
      }
    } else {

      auto subsurf_mesh = S.GetMesh(domain_);
      for (int c=0; c!=res_c.MyLength(); ++c) {
        double sum = 0;
        for (auto i : subsurf_mesh->cells_of_column(c)) {
          sum += dep_c[0][i] * cv[0][i];
        }
        res_c[0][c] = coef_ * sum / surf_cv[0][c];
      }
    }

  } else {
    if (molar_dens_key_ != "") {
      const Epetra_MultiVector& dens = *S.Get<CompositeVector>(molar_dens_key_, tag)
        .ViewComponent("cell",false);

      auto subsurf_mesh = S.GetMesh(domain_);
      for (int c=0; c!=res_c.MyLength(); ++c) {
        double sum = 0;
        for (auto i : subsurf_mesh->cells_of_column(c)) {
          sum += dep_c[0][i] / dens[0][i];
        }
        res_c[0][c] = coef_ * sum;
      }
    } else {

      auto subsurf_mesh = S.GetMesh(domain_);
      for (int c=0; c!=res_c.MyLength(); ++c) {
        double sum = 0;
        for (auto i : subsurf_mesh->cells_of_column(c)) {
          sum += dep_c[0][i];
        }
        res_c[0][c] = coef_ * sum;
      }
    }
  }
}


void
ColumnSumEvaluator::EvaluatePartialDerivative_(const State& S,
               const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(false);
}


// Custom EnsureCompatibility forces this to be updated once.
bool
ColumnSumEvaluator::Update(State& S, const Key& request)
{
  bool changed = EvaluatorSecondaryMonotypeCV::Update(S,request);

  if (!updated_once_) {
    Update_(S);
    updated_once_ = true;
    return true;
  }
  return changed;
}

void
ColumnSumEvaluator::EnsureCompatibility(State& S)
{
  Tag tag = my_keys_.front().second;
  Key key = my_keys_.front().first;

  auto& my_fac = S.Require<CompositeVector,CompositeVectorSpace>(key, tag, key);

  if (my_fac.Mesh() != Teuchos::null) {
    // Recurse into the tree to propagate info to leaves.
    for (const auto& dep : dependencies_) {
      S.RequireEvaluator(dep.first, dep.second).EnsureCompatibility(S);
    }
  }

  // check plist for vis or checkpointing control
  EvaluatorSecondaryMonotypeCV::EnsureCompatibility_Flags_(S);
}


} //namespace
} //namespace
