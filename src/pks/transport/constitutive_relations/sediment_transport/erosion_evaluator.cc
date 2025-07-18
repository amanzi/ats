/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*

 Determining erosion rate into sediment transport
  
*/

#include "erosion_evaluator.hh"

namespace Amanzi {

ErosionRateEvaluator ::ErosionRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);

  velocity_key_ = Keys::readKey(plist_, domain_name, "velocity", "velocity");
  Key pres_key = Keys::readKey(plist_, domain_name, "pressure", "pressure");

  tau_e_ = plist_.get<double>("critical shear stress");
  Qe_0_ = plist_.get<double>("empirical coefficient");
  gamma_ = plist_.get<double>("specific weight of water");
  Cf_ = plist_.get<double>("drag coefficient");

  dependencies_.insert(KeyTag{ velocity_key_, tag});
  dependencies_.insert(KeyTag{ pres_key, Tags::NEXT });
}


ErosionRateEvaluator ::ErosionRateEvaluator(const ErosionRateEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other), velocity_key_(other.velocity_key_)
{
  tau_e_ = other.tau_e_;
  Qe_0_ = other.Qe_0_;
  gamma_ = other.gamma_;
  Cf_ = other.Cf_;
  // lambda_ = other.lambda_;
  // umax_ = other.umax_;
  // xi_ = other.xi_;

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

  double max_tau = 0.0;
  double max_v2 = 0.0;

  for (int c = 0; c < vel.MyLength(); c++) {
    //double tau_0 = gamma_ * lambda_ * (sqrt(vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]));
    double v2 = vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c];
    double tau_0 = gamma_ * Cf_ * v2;

    max_tau = std::max(tau_0, max_tau);
    max_v2 = std::max(v2, max_v2);

    if (tau_0 > tau_e_) {
      result_c[0][c] = Qe_0_ * (tau_0 / tau_e_ - 1);
    } else {
      result_c[0][c] = 0.;
    }
  }

  std::cout<<"max_tau "<<max_tau<<" max v2 "<<max_v2<<"\n";
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
