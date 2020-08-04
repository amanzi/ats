/*
  The richards energy evaluator is an algebraic evaluator of a given model.
Richards energy: the standard form as a function of liquid saturation and specific internal energy.  
  Generated via evaluator_generator.
*/

#include "richards_energy_evaluator.hh"
#include "richards_energy_model.hh"

namespace Amanzi {
namespace Energy {
namespace Relations {

// Constructor from ParameterList
RichardsEnergyEvaluator::RichardsEnergyEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("richards_energy parameters");
  model_ = Teuchos::rcp(new RichardsEnergyModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
RichardsEnergyEvaluator::RichardsEnergyEvaluator(const RichardsEnergyEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    phi_key_(other.phi_key_),
    phi0_key_(other.phi0_key_),
    sl_key_(other.sl_key_),
    nl_key_(other.nl_key_),
    ul_key_(other.ul_key_),
    rho_r_key_(other.rho_r_key_),
    ur_key_(other.ur_key_),
    cv_key_(other.cv_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
RichardsEnergyEvaluator::Clone() const
{
  return Teuchos::rcp(new RichardsEnergyEvaluator(*this));
}


// Initialize by setting up dependencies
void
RichardsEnergyEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomain(my_key_);

  // - pull Keys from plist
  // dependency: porosity
  phi_key_ = plist_.get<std::string>("porosity key",
          domain_name+"porosity");
  dependencies_.insert(phi_key_);

  // dependency: base_porosity
  phi0_key_ = plist_.get<std::string>("base porosity key",
          domain_name+"base_porosity");
  dependencies_.insert(phi0_key_);

  // dependency: saturation_liquid
  sl_key_ = plist_.get<std::string>("saturation liquid key",
          domain_name+"saturation_liquid");
  dependencies_.insert(sl_key_);

  // dependency: molar_density_liquid
  nl_key_ = plist_.get<std::string>("molar density liquid key",
          domain_name+"molar_density_liquid");
  dependencies_.insert(nl_key_);

  // dependency: internal_energy_liquid
  ul_key_ = plist_.get<std::string>("internal energy liquid key",
          domain_name+"internal_energy_liquid");
  dependencies_.insert(ul_key_);

  // dependency: density_rock
  rho_r_key_ = plist_.get<std::string>("density rock key",
          domain_name+"density_rock");
  dependencies_.insert(rho_r_key_);

  // dependency: internal_energy_rock
  ur_key_ = plist_.get<std::string>("internal energy rock key",
          domain_name+"internal_energy_rock");
  dependencies_.insert(ur_key_);

  // dependency: cell_volume
  cv_key_ = plist_.get<std::string>("cell volume key",
          domain_name+"cell_volume");
  dependencies_.insert(cv_key_);
}


void
RichardsEnergyEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> phi = S->GetFieldData(phi_key_);
Teuchos::RCP<const CompositeVector> phi0 = S->GetFieldData(phi0_key_);
Teuchos::RCP<const CompositeVector> sl = S->GetFieldData(sl_key_);
Teuchos::RCP<const CompositeVector> nl = S->GetFieldData(nl_key_);
Teuchos::RCP<const CompositeVector> ul = S->GetFieldData(ul_key_);
Teuchos::RCP<const CompositeVector> rho_r = S->GetFieldData(rho_r_key_);
Teuchos::RCP<const CompositeVector> ur = S->GetFieldData(ur_key_);
Teuchos::RCP<const CompositeVector> cv = S->GetFieldData(cv_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
    const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
    const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
    const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
    const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
    const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->Energy(phi_v[0][i], phi0_v[0][i], sl_v[0][i], nl_v[0][i], ul_v[0][i], rho_r_v[0][i], ur_v[0][i], cv_v[0][i]);
    }
  }
}


void
RichardsEnergyEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> phi = S->GetFieldData(phi_key_);
Teuchos::RCP<const CompositeVector> phi0 = S->GetFieldData(phi0_key_);
Teuchos::RCP<const CompositeVector> sl = S->GetFieldData(sl_key_);
Teuchos::RCP<const CompositeVector> nl = S->GetFieldData(nl_key_);
Teuchos::RCP<const CompositeVector> ul = S->GetFieldData(ul_key_);
Teuchos::RCP<const CompositeVector> rho_r = S->GetFieldData(rho_r_key_);
Teuchos::RCP<const CompositeVector> ur = S->GetFieldData(ur_key_);
Teuchos::RCP<const CompositeVector> cv = S->GetFieldData(cv_key_);

  if (wrt_key == phi_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDPorosity(phi_v[0][i], phi0_v[0][i], sl_v[0][i], nl_v[0][i], ul_v[0][i], rho_r_v[0][i], ur_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == phi0_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDBasePorosity(phi_v[0][i], phi0_v[0][i], sl_v[0][i], nl_v[0][i], ul_v[0][i], rho_r_v[0][i], ur_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == sl_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDSaturationLiquid(phi_v[0][i], phi0_v[0][i], sl_v[0][i], nl_v[0][i], ul_v[0][i], rho_r_v[0][i], ur_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == nl_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDMolarDensityLiquid(phi_v[0][i], phi0_v[0][i], sl_v[0][i], nl_v[0][i], ul_v[0][i], rho_r_v[0][i], ur_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == ul_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDInternalEnergyLiquid(phi_v[0][i], phi0_v[0][i], sl_v[0][i], nl_v[0][i], ul_v[0][i], rho_r_v[0][i], ur_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == rho_r_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDDensityRock(phi_v[0][i], phi0_v[0][i], sl_v[0][i], nl_v[0][i], ul_v[0][i], rho_r_v[0][i], ur_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == ur_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDInternalEnergyRock(phi_v[0][i], phi0_v[0][i], sl_v[0][i], nl_v[0][i], ul_v[0][i], rho_r_v[0][i], ur_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == cv_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDCellVolume(phi_v[0][i], phi0_v[0][i], sl_v[0][i], nl_v[0][i], ul_v[0][i], rho_r_v[0][i], ur_v[0][i], cv_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
