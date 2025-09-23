/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

*/

#include "iem_evaluator.hh"
#include "iem_factory.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {


IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  AMANZI_ASSERT(plist_.isSublist("IEM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("IEM parameters");
  IEMFactory fac;
  iem_ = fac.createIEM(sublist);

  InitializeFromPlist_();
}


IEMEvaluator::IEMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEM>& iem)
  : EvaluatorSecondaryMonotypeCV(plist), iem_(iem)
{
  InitializeFromPlist_();
}


Teuchos::RCP<Evaluator>
IEMEvaluator::Clone() const
{
  return Teuchos::rcp(new IEMEvaluator(*this));
}


void
IEMEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- temperature
  temp_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });
}


void
IEMEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = iem_->InternalEnergy(temp_v[0][i]);
    }
  }
}


void
IEMEvaluator::EvaluatePartialDerivative_(const State& S,
                                         const Key& wrt_key,
                                         const Tag& wrt_tag,
                                         const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(wrt_key == temp_key_);
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = iem_->DInternalEnergyDT(temp_v[0][i]);
    }
  }
}


} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
