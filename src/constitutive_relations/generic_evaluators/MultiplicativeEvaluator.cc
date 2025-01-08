/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

#include "MultiplicativeEvaluator.hh"

namespace Amanzi {
namespace Relations {

MultiplicativeEvaluator::MultiplicativeEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  if (dependencies_.size() == 0) {
    Errors::Message message;
    message << "MultiplicativeEvaluator: for " << my_keys_[0].first
            << " was provided no dependencies";
    throw(message);
  }

  // deals with Multiplying e.g. area fractions.  Note that 90% of
  // multiplicative evaluator usages will NOT provide dof info, and therefore
  // will be the same structure as my_key.  Not having this info is stronger,
  // because we can set the stencil.  Prefer not to provide it.
  any_dof_provided_ = false;
  for (const auto& dep : dependencies_) {
    if (plist_.isParameter(dep.first + " degree of freedom")) {
      dofs_.push_back(plist_.get<int>(dep.first + " degree of freedom"));
      dof_provided_.push_back(true);
      any_dof_provided_ = true;
    } else {
      dofs_.push_back(0);
      dof_provided_.push_back(false);
    }
  }

  coef_ = plist_.get<double>("coefficient", 1.0);
  positive_ = plist_.get<bool>("enforce positivity", false);
}


Teuchos::RCP<Evaluator>
MultiplicativeEvaluator::Clone() const
{
  return Teuchos::rcp(new MultiplicativeEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
MultiplicativeEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(dependencies_.size() >= 1);
  result[0]->PutScalar(coef_);

  for (const auto& lcv_name : *result[0]) {
    // note, this multiply is done with Vectors, not MultiVectors, to allow DoFs
    auto& res_c = *(result[0]->ViewComponent(lcv_name, false));
    int i = 0;
    for (const auto& key_tag : dependencies_) {
      const auto& dep_v =
        *(*S.Get<CompositeVector>(key_tag.first, key_tag.second).ViewComponent(lcv_name, false))(
          dofs_[i]);
      res_c.Multiply(1, dep_v, res_c, 0.);
      i++;
    }

    if (positive_) {
      for (int c = 0; c != res_c.MyLength(); ++c) {
        for (int i = 0; i != res_c.NumVectors(); ++i) { res_c[i][c] = std::max(res_c[i][c], 0.); }
      }
    }
  }
}

void
MultiplicativeEvaluator::EvaluatePartialDerivative_(const State& S,
                                                    const Key& wrt_key,
                                                    const Tag& wrt_tag,
                                                    const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(coef_);

  for (const auto& lcv_name : *result[0]) {
    // note, this multiply is done with Vectors, not MultiVectors, to allow DoFs
    auto& res_c = *(result[0]->ViewComponent(lcv_name, false));
    int i = 0;
    for (const auto& key_tag : dependencies_) {
      if ((key_tag.first != wrt_key) || (key_tag.second != wrt_tag)) {
        const auto& dep_v =
          *(*S.Get<CompositeVector>(key_tag.first, key_tag.second).ViewComponent(lcv_name, false))(
            dofs_[i]);
        res_c.Multiply(1, dep_v, res_c, 0.);
        i++;
      }
    }

    if (positive_) {
      const auto& value_c = *S.Get<CompositeVector>(my_keys_.front().first, my_keys_.front().second)
                               .ViewComponent(lcv_name, false);
      for (int c = 0; c != res_c.MyLength(); ++c) {
        for (int i = 0; i != res_c.NumVectors(); ++i) {
          if (value_c[i][c] == 0) { res_c[i][c] = 0.; }
        }
      }
    }
  }
}


void
MultiplicativeEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  if (any_dof_provided_) {
    Key my_key = my_keys_.front().first;
    Tag my_tag = my_keys_.front().second;

    // If my requirements have not yet been set, we'll have to hope they
    // get set by someone later.  For now just defer.
    auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_key, my_tag);
    if (my_fac.Mesh() != Teuchos::null) {
      // Create an unowned factory to check my dependencies.
      CompositeVectorSpace dep_fac(my_fac);
      dep_fac.SetOwned(false);

      // Loop over my dependencies, ensuring they meet the requirements.
      int i = 0;
      for (const auto& key_tag : dependencies_) {
        if (key_tag.first == my_key && key_tag.second == my_tag) {
          Errors::Message msg;
          msg << "Evaluator for key \"" << my_key << " @ " << my_tag.get()
              << "\" depends upon itself.";
          Exceptions::amanzi_throw(msg);
        }
        auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second);

        if (!dof_provided_[i]) {
          // this must be a 1-dof entry of the same shape as my_key
          fac.Update(dep_fac);
        } else {
          // may have different shape, just update mesh
          // just have to hope that # of dofs and stencil is set later
          fac.SetMesh(dep_fac.Mesh());
        }
      }
    }

  } else {
    // the standard ensure is better
    return EvaluatorSecondaryMonotypeCV::EnsureCompatibility_ToDeps_(S);
  }
}

} // namespace Relations
} // namespace Amanzi
