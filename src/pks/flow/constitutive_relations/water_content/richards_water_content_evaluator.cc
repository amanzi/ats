/*
  The richards water content evaluator is an algebraic evaluator of a given model.
Richards water content evaluator: the standard form as a function of liquid saturation.  
  Generated via evaluator_generator.
*/

#include "richards_water_content_evaluator.hh"
#include "richards_water_content_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
RichardsWaterContentEvaluator::RichardsWaterContentEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("richards_water_content parameters");
  model_ = Teuchos::rcp(new RichardsWaterContentModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
RichardsWaterContentEvaluator::RichardsWaterContentEvaluator(const RichardsWaterContentEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    phi_key_(other.phi_key_),
    sl_key_(other.sl_key_),
    nl_key_(other.nl_key_),
    cv_key_(other.cv_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
RichardsWaterContentEvaluator::Clone() const
{
  return Teuchos::rcp(new RichardsWaterContentEvaluator(*this));
}


// Initialize by setting up dependencies
void
RichardsWaterContentEvaluator::InitializeFromPlist_()
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

  // dependency: cell_volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(cv_key_);
}


void
RichardsWaterContentEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> phi = S->GetFieldData(phi_key_);
Teuchos::RCP<const CompositeVector> sl = S->GetFieldData(sl_key_);
Teuchos::RCP<const CompositeVector> nl = S->GetFieldData(nl_key_);
Teuchos::RCP<const CompositeVector> cv = S->GetFieldData(cv_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
    const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->WaterContent(phi_v[0][i], sl_v[0][i], nl_v[0][i], cv_v[0][i]);
    }
  }
}


void
RichardsWaterContentEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> phi = S->GetFieldData(phi_key_);
Teuchos::RCP<const CompositeVector> sl = S->GetFieldData(sl_key_);
Teuchos::RCP<const CompositeVector> nl = S->GetFieldData(nl_key_);
Teuchos::RCP<const CompositeVector> cv = S->GetFieldData(cv_key_);

  if (wrt_key == phi_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDPorosity(phi_v[0][i], sl_v[0][i], nl_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == sl_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDSaturationLiquid(phi_v[0][i], sl_v[0][i], nl_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == nl_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDMolarDensityLiquid(phi_v[0][i], sl_v[0][i], nl_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == cv_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDCellVolume(phi_v[0][i], sl_v[0][i], nl_v[0][i], cv_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
