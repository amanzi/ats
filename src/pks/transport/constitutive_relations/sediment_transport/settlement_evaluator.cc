/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Determining rate of settlement of sediment
*/

#include "settlement_evaluator.hh"

namespace Amanzi {

SettlementRateEvaluator::SettlementRateEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain_name = Keys::getDomain(my_keys_.front().first);

  velocity_key_ = Keys::readKey(plist_, domain_name, "velocity", "velocity");
  sediment_key_ = Keys::readKey(plist_, domain_name, "sediment", "sediment");

  dependencies_.insert(KeyTag{ sediment_key_, tag });
  dependencies_.insert(KeyTag{ velocity_key_, tag });

  tau_d_ = plist_.get<double>("critical shear stress");
  ws_ = plist_.get<double>("settling velocity");
  gamma_ = plist_.get<double>("specific weight of water");

  Cf_ = plist_.get<double>("drag coefficient");
}


Teuchos::RCP<Evaluator>
SettlementRateEvaluator::Clone() const
{
  return Teuchos::rcp(new SettlementRateEvaluator(*this));
}


void
SettlementRateEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& vel =
    *S.GetPtr<CompositeVector>(velocity_key_, tag)->ViewComponent("cell");
  const Epetra_MultiVector& tcc = *S.Get<CompositeVector>(sediment_key_, tag).ViewComponent("cell");
  Epetra_MultiVector& result_c = *result[0]->ViewComponent("cell");

  sediment_density_ = S.Get<double>("sediment_density", tag);

  for (int c = 0; c < result_c.MyLength(); c++) {
    double tau_0 = gamma_ * Cf_ * (vel[0][c] * vel[0][c] + vel[1][c] * vel[1][c]);

    if (tau_0 < tau_d_) {
      result_c[0][c] = sediment_density_ * ws_ * std::min(tcc[0][c], 0.5) * (1 - tau_0 / tau_d_);
    } else {
      result_c[0][c] = 0.;
    }
  }
}

void
SettlementRateEvaluator::EvaluatePartialDerivative_(const State& S,
                                                    const Key& wrt_key,
                                                    const Tag& wrt_tag,
                                                    const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(0);
}

} // namespace Amanzi
