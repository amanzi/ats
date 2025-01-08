/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  PCLiquidEvaluator is the interface between state/data and the model, a PC relation.

*/

#include "pc_liq_atm.hh"
#include "pc_liquid_evaluator.hh"

namespace Amanzi {
namespace Flow {

PCLiquidEvaluator::PCLiquidEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  // dependencies
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- pressure
  pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(KeyTag{ pres_key_, tag });

  // Construct my PCLiquid model
  model_ = Teuchos::rcp(new PCLiqAtm(plist_.sublist("capillary pressure model parameters")));
};


Teuchos::RCP<Evaluator>
PCLiquidEvaluator::Clone() const
{
  return Teuchos::rcp(new PCLiquidEvaluator(*this));
}


void
PCLiquidEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);
  const double& p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);

  // evaluate pc
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int id = 0; id != count; ++id) {
      result_v[0][id] = model_->CapillaryPressure(pres_v[0][id], p_atm);
    }
  }
}


void
PCLiquidEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(wrt_key == pres_key_);

  // Pull dependencies out of state.
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);
  const double& p_atm = S.Get<double>("atmospheric_pressure", Tags::DEFAULT);

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& pres_v = *(pres->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int id = 0; id != count; ++id) {
      result_v[0][id] = model_->DCapillaryPressureDp(pres_v[0][id], p_atm);
    }
  }
}

} // namespace Flow
} // namespace Amanzi
