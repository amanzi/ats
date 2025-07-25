/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Determining rate for trapping sediment
*/

#include "trapping_evaluator.hh"

namespace Amanzi {

TrappingRateEvaluator ::TrappingRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);

  velocity_key_ = Keys::readKey(plist_, domain_name, "velocity", "velocity");
  dependencies_.insert(KeyTag{ velocity_key_, tag });

  sediment_key_ = Keys::readKey(plist_, domain_name, "sediment", "sediment");
  dependencies_.insert(KeyTag{ sediment_key_, tag });

  ponded_depth_key_ = Keys::readKey(plist_, domain_name, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ ponded_depth_key_, tag });

  stem_diameter_key_ = Keys::readKey(plist_, domain_name, "stem_diameter", "stem_diameter");
  dependencies_.insert(KeyTag{ stem_diameter_key_, tag });

  stem_height_key_ = Keys::readKey(plist_, domain_name, "stem_height", "stem_height");
  dependencies_.insert(KeyTag{ stem_height_key_, tag });

  stem_density_key_ = Keys::readKey(plist_, domain_name, "stem_density", "stem_density");
  dependencies_.insert(KeyTag{ stem_density_key_, tag });

  // Please put units on all of these! --ETC
  visc_ = plist_.get<double>("kinematic viscosity");
  d_p_ = plist_.get<double>("particle diameter");
  alpha_ = plist_.get<double>("alpha");
  beta_ = plist_.get<double>("beta");
  gamma_ = plist_.get<double>("gamma");
}


Teuchos::RCP<Evaluator>
TrappingRateEvaluator::Clone() const
{
  return Teuchos::rcp(new TrappingRateEvaluator(*this));
}


void
TrappingRateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  sediment_density_ = S.Get<double>("sediment_density", tag);
  const Epetra_MultiVector& vel = *S.Get<CompositeVector>(velocity_key_, tag).ViewComponent("cell");
  const Epetra_MultiVector& tcc = *S.Get<CompositeVector>(sediment_key_, tag).ViewComponent("cell");
  const Epetra_MultiVector& depth =
    *S.Get<CompositeVector>(ponded_depth_key_, tag).ViewComponent("cell");
  const Epetra_MultiVector& bio_n =
    *S.Get<CompositeVector>(stem_density_key_, tag).ViewComponent("cell");
  const Epetra_MultiVector& bio_d =
    *S.Get<CompositeVector>(stem_diameter_key_, tag).ViewComponent("cell");
  const Epetra_MultiVector& bio_h =
    *S.Get<CompositeVector>(stem_height_key_, tag).ViewComponent("cell");
  Epetra_MultiVector& result_c = *result[0]->ViewComponent("cell");

  result_c.PutScalar(0.);

  for (int c = 0; c < result_c.MyLength(); c++) {
    for (int j = 0; j < bio_n.NumVectors(); j++) {
      double h_s = bio_h[j][c];
      double d_s = bio_d[j][c];
      double n_s = bio_n[j][c];

      double u_abs = sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]);
      double eps;
      if (d_s > 1e-12) {
        eps = alpha_ * std::pow(u_abs * d_s / visc_, beta_) * std::pow(d_p_ / d_s, gamma_);
      } else {
        eps = 0.;
      }

      result_c[0][c] +=
        sediment_density_ * tcc[0][c] * u_abs * eps * d_s * n_s * std::min(depth[0][c], h_s);
    }
  }
}


void
TrappingRateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                  const Key& wrt_key,
                                                  const Tag& wrt_tag,
                                                  const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
}

} // namespace Amanzi
