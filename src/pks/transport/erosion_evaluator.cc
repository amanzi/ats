/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*


*/

#include "erosion_evaluator.hh"
#include "boost/math/constants/constants.hpp"

namespace Amanzi {
  
ErosionRateEvaluator ::ErosionRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);

  velocity_key_ = plist_.get<std::string>("velocity key", Keys::getKey(domain_name, "velocity"));

  tau_e_ = plist_.get<double>("critical shear stress");
  Qe_0_ = plist_.get<double>("empirical coefficient");
  gamma_ = plist_.get<double>("specific weight of water");
  umax_ = plist_.get<double>("max current");
  xi_ = plist_.get<double>("Chezy parameter");
  Cf_ = plist_.get<double>("drag coefficient");

  double pi = boost::math::constants::pi<double>();

  lambda_ = 8. / (3 * pi) * (umax_ / (xi_ * xi_));

  dependencies_.insert(KeyTag{ "surface-pressure", tag });
}


ErosionRateEvaluator ::ErosionRateEvaluator(const ErosionRateEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other), velocity_key_(other.velocity_key_)
{
  tau_e_ = other.tau_e_;
  Qe_0_ = other.Qe_0_;
  gamma_ = other.gamma_;
  lambda_ = other.lambda_;
  umax_ = other.umax_;
  xi_ = other.xi_;
  Cf_ = other.Cf_;
}


Teuchos::RCP<Evaluator>
ErosionRateEvaluator ::Clone() const
{
  return Teuchos::rcp(new ErosionRateEvaluator(*this));
}


void
ErosionRateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& vel =
    *S.GetPtr<CompositeVector>(velocity_key_, tag)->ViewComponent("cell");
  Epetra_MultiVector& result_c = *result[0]->ViewComponent("cell");

  for (int c = 0; c < vel.MyLength(); c++) {
    //double tau_0 = gamma_ * lambda_ * (sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]));
    double tau_0 = gamma_ * Cf_ *
                   (sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]) *
                    sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]));
    if (tau_0 > tau_e_) {
      result_c[0][c] = Qe_0_ * (tau_0 / tau_e_ - 1);
    } else {
      result_c[0][c] = 0.;
    }
  }
}

void
ErosionRateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                 const Key& wrt_key,
                                                 const Tag& wrt_tag,
                                                 const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
}

} // namespace Amanzi
