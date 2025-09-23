/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluator for determining height( rho, head )

*/

#include "effective_height_model.hh"
#include "effective_height_evaluator.hh"


namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

EffectiveHeightEvaluator::EffectiveHeightEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  auto domain_name = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  // my dependencies
  height_key_ = Keys::readKey(plist_, domain_name, "height key", "ponded_depth");
  dependencies_.insert(KeyTag{ height_key_, tag });

  // model
  Teuchos::ParameterList model_plist = plist_.sublist("effective height model parameters");
  model_ = Teuchos::rcp(new EffectiveHeightModel(model_plist));
}


Teuchos::RCP<Evaluator>
EffectiveHeightEvaluator::Clone() const
{
  return Teuchos::rcp(new EffectiveHeightEvaluator(*this));
}

void
EffectiveHeightEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> height = S.GetPtr<CompositeVector>(height_key_, tag);

  // evaluate p_s / p_atm
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& height_v = *(height->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int id = 0; id != count; ++id) {
      result_v[0][id] = model_->EffectiveHeight(height_v[0][id]);
    }
  }
}


void
EffectiveHeightEvaluator::EvaluatePartialDerivative_(const State& S,
                                                     const Key& wrt_key,
                                                     const Tag& wrt_tag,
                                                     const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(wrt_key == height_key_);

  // Pull dependencies out of state.
  auto tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> height = S.GetPtr<CompositeVector>(height_key_, tag);

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& height_v = *(height->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int id = 0; id != count; ++id) {
      result_v[0][id] = model_->DEffectiveHeightDHeight(height_v[0][id]);
    }
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
