/*
  The interfrost denergy_dtemperature evaluator is an algebraic evaluator of a given model.
Interfrost water content portion sl.  
  Generated via evaluator_generator.
*/

#include "interfrost_denergy_dtemperature_evaluator.hh"
#include "interfrost_denergy_dtemperature_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
InterfrostDenergyDtemperatureEvaluator::InterfrostDenergyDtemperatureEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("interfrost_denergy_dtemperature parameters");
  model_ = Teuchos::rcp(new InterfrostDenergyDtemperatureModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
InterfrostDenergyDtemperatureEvaluator::InterfrostDenergyDtemperatureEvaluator(const InterfrostDenergyDtemperatureEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    phi_key_(other.phi_key_),
    sl_key_(other.sl_key_),
    nl_key_(other.nl_key_),
    si_key_(other.si_key_),
    ni_key_(other.ni_key_),
    rhos_key_(other.rhos_key_),
    T_key_(other.T_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
InterfrostDenergyDtemperatureEvaluator::Clone() const
{
  return Teuchos::rcp(new InterfrostDenergyDtemperatureEvaluator(*this));
}


// Initialize by setting up dependencies
void
InterfrostDenergyDtemperatureEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomain(my_key_);

  // - pull Keys from plist
  // dependency: porosity
  phi_key_ = Keys::readKey(plist_, domain_name, "porosity", "porosity");
  dependencies_.insert(phi_key_);

  // dependency: saturation_liquid
  sl_key_ = Keys::readKey(plist_, domain_name, "saturation liquid", "saturation_liquid");
  dependencies_.insert(sl_key_);

  // dependency: molar_density_liquid
  nl_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(nl_key_);

  // dependency: saturation_ice
  si_key_ = Keys::readKey(plist_, domain_name, "saturation ice", "saturation_ice");
  dependencies_.insert(si_key_);

  // dependency: molar_density_ice
  ni_key_ = Keys::readKey(plist_, domain_name, "molar density ice", "molar_density_ice");
  dependencies_.insert(ni_key_);

  // dependency: density_rock
  rhos_key_ = Keys::readKey(plist_, domain_name, "density rock", "density_rock");
  dependencies_.insert(rhos_key_);

  // dependency: temperature
  T_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(T_key_);
}


void
InterfrostDenergyDtemperatureEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> phi = S->GetFieldData(phi_key_);
Teuchos::RCP<const CompositeVector> sl = S->GetFieldData(sl_key_);
Teuchos::RCP<const CompositeVector> nl = S->GetFieldData(nl_key_);
Teuchos::RCP<const CompositeVector> si = S->GetFieldData(si_key_);
Teuchos::RCP<const CompositeVector> ni = S->GetFieldData(ni_key_);
Teuchos::RCP<const CompositeVector> rhos = S->GetFieldData(rhos_key_);
Teuchos::RCP<const CompositeVector> T = S->GetFieldData(T_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
    const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
    const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
    const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
    const Epetra_MultiVector& rhos_v = *rhos->ViewComponent(*comp, false);
    const Epetra_MultiVector& T_v = *T->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->DEnergyDTCoef(phi_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], rhos_v[0][i], T_v[0][i]);
    }
  }
}


void
InterfrostDenergyDtemperatureEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> phi = S->GetFieldData(phi_key_);
Teuchos::RCP<const CompositeVector> sl = S->GetFieldData(sl_key_);
Teuchos::RCP<const CompositeVector> nl = S->GetFieldData(nl_key_);
Teuchos::RCP<const CompositeVector> si = S->GetFieldData(si_key_);
Teuchos::RCP<const CompositeVector> ni = S->GetFieldData(ni_key_);
Teuchos::RCP<const CompositeVector> rhos = S->GetFieldData(rhos_key_);
Teuchos::RCP<const CompositeVector> T = S->GetFieldData(T_key_);

  if (wrt_key == phi_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& rhos_v = *rhos->ViewComponent(*comp, false);
      const Epetra_MultiVector& T_v = *T->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DDEnergyDTCoefDPorosity(phi_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], rhos_v[0][i], T_v[0][i]);
      }
    }

  } else if (wrt_key == sl_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& rhos_v = *rhos->ViewComponent(*comp, false);
      const Epetra_MultiVector& T_v = *T->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DDEnergyDTCoefDSaturationLiquid(phi_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], rhos_v[0][i], T_v[0][i]);
      }
    }

  } else if (wrt_key == nl_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& rhos_v = *rhos->ViewComponent(*comp, false);
      const Epetra_MultiVector& T_v = *T->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DDEnergyDTCoefDMolarDensityLiquid(phi_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], rhos_v[0][i], T_v[0][i]);
      }
    }

  } else if (wrt_key == si_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& rhos_v = *rhos->ViewComponent(*comp, false);
      const Epetra_MultiVector& T_v = *T->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DDEnergyDTCoefDSaturationIce(phi_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], rhos_v[0][i], T_v[0][i]);
      }
    }

  } else if (wrt_key == ni_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& rhos_v = *rhos->ViewComponent(*comp, false);
      const Epetra_MultiVector& T_v = *T->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DDEnergyDTCoefDMolarDensityIce(phi_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], rhos_v[0][i], T_v[0][i]);
      }
    }

  } else if (wrt_key == rhos_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& rhos_v = *rhos->ViewComponent(*comp, false);
      const Epetra_MultiVector& T_v = *T->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DDEnergyDTCoefDDensityRock(phi_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], rhos_v[0][i], T_v[0][i]);
      }
    }

  } else if (wrt_key == T_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& rhos_v = *rhos->ViewComponent(*comp, false);
      const Epetra_MultiVector& T_v = *T->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DDEnergyDTCoefDTemperature(phi_v[0][i], sl_v[0][i], nl_v[0][i], si_v[0][i], ni_v[0][i], rhos_v[0][i], T_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
