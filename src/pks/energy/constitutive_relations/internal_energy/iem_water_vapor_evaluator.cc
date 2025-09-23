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

#include "iem_water_vapor_evaluator.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {

IEMWaterVaporEvaluator::IEMWaterVaporEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  // defaults work fine, this sublist need not exist
  Teuchos::ParameterList sublist = plist.sublist("IEM parameters");
  iem_ = Teuchos::rcp(new IEMWaterVapor(sublist));
  InitializeFromPlist_();
}


IEMWaterVaporEvaluator::IEMWaterVaporEvaluator(Teuchos::ParameterList& plist,
                                               const Teuchos::RCP<IEMWaterVapor>& iem)
  : EvaluatorSecondaryMonotypeCV(plist), iem_(iem)
{
  InitializeFromPlist_();
}


Teuchos::RCP<Evaluator>
IEMWaterVaporEvaluator::Clone() const
{
  return Teuchos::rcp(new IEMWaterVaporEvaluator(*this));
}


void
IEMWaterVaporEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies.
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- temperature
  temp_key_ = Keys::readKey(plist_, domain_name, "temperature key", "temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });

  // -- molar fraction of water vapor in the gaseous phase
  mol_frac_key_ = Keys::readKey(plist_, domain_name, "vapor molar fraction key", "mol_frac_gas");
  dependencies_.insert(KeyTag{ mol_frac_key_, tag });
}


void
IEMWaterVaporEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);
  Teuchos::RCP<const CompositeVector> mol_frac = S.GetPtr<CompositeVector>(mol_frac_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = iem_->InternalEnergy(temp_v[0][i], molfrac_v[0][i]);
    }
  }
}


void
IEMWaterVaporEvaluator::EvaluatePartialDerivative_(const State& S,
                                                   const Key& wrt_key,
                                                   const Tag& wrt_tag,
                                                   const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);
  Teuchos::RCP<const CompositeVector> mol_frac = S.GetPtr<CompositeVector>(mol_frac_key_, tag);

  if (wrt_key == temp_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
      const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = iem_->DInternalEnergyDT(temp_v[0][i], molfrac_v[0][i]);
      }
    }
  } else if (wrt_key == mol_frac_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
      const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = iem_->DInternalEnergyDomega(temp_v[0][i], molfrac_v[0][i]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}


} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
