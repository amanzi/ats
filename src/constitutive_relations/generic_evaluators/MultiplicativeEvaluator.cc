/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "MultiplicativeEvaluator.hh"

namespace Amanzi {
namespace Relations {

MultiplicativeEvaluator::MultiplicativeEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
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
    if (plist_.isParameter(dep+" degree of freedom")) {
      dofs_.push_back(plist_.get<int>(dep+" degree of freedom"));
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
MultiplicativeEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(dependencies_.size() > 1);
  result[0]->PutScalar(coef_);

  for (const auto& lcv_name : *result[0]) {
    // note, this multiply is done with Vectors, not MultiVectors, to allow DoFs
    auto& res_c = *(*result[0]->ViewComponent(lcv_name, false))(0);
    int i = 0;
    for (const auto& key_tag : dependencies_) {
      const auto& dep_v = *(*S->Get<CompositeVector>(key_tag.first, key_tag.second)
                           .ViewComponent(lcv_name, false))(dofs_[i]);
      res_c.Multiply(1, res_c, dep_v, 0.);
      i++;
    }

    if (positive_) {
      for (int c=0; c!=res_c.MyLength(); ++c) {
        res_c[c] = std::max(res_c[c], 0.);
      }
    }
  }
}

void
MultiplicativeEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag,
        const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(coef_);

  for (const auto& lcv_name : *result[0]) {
    // note, this multiply is done with Vectors, not MultiVectors, to allow DoFs
    auto& res_c = *(*result[0]->ViewComponent(lcv_name, false))(0);
    int i = 0;
    for (const auto& key_tag : dependencies_) {
      if (key != wrt_key) {
        const auto& dep_v = *(*S->Get<CompositeVector>(key_tag.first, key_tag.second)
                .ViewComponent(lcv_name, false))(dofs_[i]);
        res_c.Multiply(1, res_c, dep_v, 0.);
        i++;
      }
    }
    if (positive_) {
      for (int c=0; c!=res_c.MyLength(); ++c) {
        res_c[c] = std::max(res_c[c], 0.);
      }
    }
  }
}


void
MultiplicativeEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S)
{
  if (any_dof_provided_) {
    // Ensure my field exists.  Requirements should be already set.
    AMANZI_ASSERT(my_key_ != std::string(""));
    Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

    // check plist for vis or checkpointing control
    bool io_my_key = plist_.get<bool>("visualize", true);
    S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>("checkpoint", false);
    S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

    // If my requirements have not yet been set, we'll have to hope they
    // get set by someone later.  For now just defer.
    if (my_fac->Mesh() != Teuchos::null) {
      // Create an unowned factory to check my dependencies.
      Teuchos::RCP<CompositeVectorSpace> dep_fac =
        Teuchos::rcp(new CompositeVectorSpace(*my_fac));
      dep_fac->SetOwned(false);

      // Loop over my dependencies, ensuring they meet the requirements.
      int i = 0;
      for (const auto& key : dependencies_) {
        if (key == my_key_) {
          Errors::Message msg;
          msg << "Evaluator for key \"" << my_key_ << "\" depends upon itself.";
          Exceptions::amanzi_throw(msg);
        }
        Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(key);

        if (!dof_provided_[i]) {
          // this must be a 1-dof entry of the same shape as my_key
          fac->Update(*dep_fac);
        } else {
          // may have different shape, just update mesh
          // just have to hope that # of dofs and stencil is set later
          fac->SetMesh(dep_fac->Mesh());
        }
      }

      // Recurse into the tree to propagate info to leaves.
      for (const auto& key : dependencies_) {
        S->RequireFieldEvaluator(key)->EnsureCompatibility(S);
      }
    }

  } else {
    // the standard ensure is better, as it checks components
    return SecondaryVariableFieldEvaluator::EnsureCompatibility(S);
  }
}


} // namespace
} // namespace

