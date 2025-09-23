/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Evaluates the conductivity of surface flow.

*/

#include "unfrozen_fraction_model.hh"
#include "unfrozen_fraction_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Flow {

UnfrozenFractionEvaluator::UnfrozenFractionEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Tag tag = my_keys_.front().second;
  Key domain = Keys::getDomain(my_keys_.front().first);

  temp_key_ = Keys::readKey(plist_, domain, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });

  // create the model, hard-coded until we have a 2nd model
  AMANZI_ASSERT(plist_.isSublist("unfrozen fraction model"));
  Teuchos::ParameterList sublist = plist_.sublist("unfrozen fraction model");
  model_ = Teuchos::rcp(new UnfrozenFractionModel(sublist));
}


Teuchos::RCP<Evaluator>
UnfrozenFractionEvaluator::Clone() const
{
  return Teuchos::rcp(new UnfrozenFractionEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
UnfrozenFractionEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = model_->UnfrozenFraction(temp_v[0][i]);
    }
  }
}


void
UnfrozenFractionEvaluator::EvaluatePartialDerivative_(const State& S,
                                                      const Key& wrt_key,
                                                      const Tag& wrt_tag,
                                                      const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  AMANZI_ASSERT(wrt_key == temp_key_);
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = model_->DUnfrozenFractionDT(temp_v[0][i]);
    }
  }
}


} // namespace Flow
} // namespace ATS_Physics
} // namespace Amanzi
