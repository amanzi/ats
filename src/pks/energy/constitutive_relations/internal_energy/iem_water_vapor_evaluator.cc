/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "iem_water_vapor_evaluator.hh"

namespace Amanzi {
namespace Energy {

IEMWaterVaporEvaluator::IEMWaterVaporEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  // defaults work fine, this sublist need not exist
  Teuchos::ParameterList sublist = plist.sublist("IEM parameters");
  iem_ = Teuchos::rcp(new IEMWaterVapor(sublist));

  InitializeFromPlist_();
}

IEMWaterVaporEvaluator::IEMWaterVaporEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<IEMWaterVapor>& iem) :
    SecondaryVariableFieldEvaluator(plist),
    iem_(iem) {
  InitializeFromPlist_();
}

IEMWaterVaporEvaluator::IEMWaterVaporEvaluator(const IEMWaterVaporEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    iem_(other.iem_),
    temp_key_(other.temp_key_),
    mol_frac_key_(other.mol_frac_key_) {}

Teuchos::RCP<FieldEvaluator>
IEMWaterVaporEvaluator::Clone() const {
  return Teuchos::rcp(new IEMWaterVaporEvaluator(*this));
}

void IEMWaterVaporEvaluator::InitializeFromPlist_() {
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("internal energy key");
  }

  // Set up my dependencies.
  Key domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temp_key_ = plist_.get<std::string>("temperature key",
                                      Keys::getKey(domain_name,std::string("temperature")));
  dependencies_.insert(temp_key_);

  // -- molar fraction of water vapor in the gaseous phase
  mol_frac_key_ = plist_.get<std::string>("vapor molar fraction key",
                                          Keys::getKey(domain_name,std::string("mol_frac_gas")));
  dependencies_.insert(mol_frac_key_);
}


void IEMWaterVaporEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);
  Teuchos::RCP<const CompositeVector> mol_frac = S->GetPtr<CompositeVector>(mol_frac_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = iem_->InternalEnergy(temp_v[0][i], molfrac_v[0][i]);
    }
  }
}


void IEMWaterVaporEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);
  Teuchos::RCP<const CompositeVector> mol_frac = S->GetPtr<CompositeVector>(mol_frac_key_);

  if (wrt_key == temp_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = iem_->DInternalEnergyDT(temp_v[0][i], molfrac_v[0][i]);
      }
    }
  } else if (wrt_key == mol_frac_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      const Epetra_MultiVector& molfrac_v = *mol_frac->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = iem_->DInternalEnergyDomega(temp_v[0][i], molfrac_v[0][i]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
