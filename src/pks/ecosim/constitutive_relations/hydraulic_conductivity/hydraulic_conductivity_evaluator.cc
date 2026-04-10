/*
  The hydraulic conductivity evaluator is an algebraic evaluator of a given model.
Richards water content evaluator: the standard form as a function of liquid saturation.
  Generated via evaluator_generator.
*/

#include "hydraulic_conductivity_evaluator.hh"
#include "hydraulic_conductivity_model.hh"

namespace Amanzi {
namespace Ecosim {
namespace Relations {

// Constructor from ParameterList
HydraulicConductivityEvaluator::HydraulicConductivityEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("hydraulic_conductivity parameters");
  model_ = Teuchos::rcp(new HydraulicConductivityModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
//Don't seem to need this
/*HydraulicConductivityEvaluator::HydraulicConductivityEvaluator(const HydraulicConductivityEvaluator& other) :
    EvaluatorSecondaryMonotypeCV(other),
    k_key_(other.k_key_),
    rho_key_(other.rho_key_),
    mu_key_(other.mu_key_),
    model_(other.model_) {}*/


// Virtual copy constructor
Teuchos::RCP<Evaluator>
HydraulicConductivityEvaluator::Clone() const
{
  return Teuchos::rcp(new HydraulicConductivityEvaluator(*this));
}


// Initialize by setting up dependencies
void
HydraulicConductivityEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  //Key domain_name = Keys::getDomain(my_key_);
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: permeability
  k_key_ = Keys::readKey(plist_, domain_name, "permeability", "permeability");
  dependencies_.insert(KeyTag{ k_key_, tag });

  // dependency: mass_density_liquid
  rho_key_ = Keys::readKey(plist_, domain_name, "mass density liquid", "mass_density_liquid");
  dependencies_.insert(KeyTag{ rho_key_, tag});

  // dependency: viscosity_liquid
  mu_key_ = Keys::readKey(plist_, domain_name, "viscosity liquid", "viscosity_liquid");
  dependencies_.insert(KeyTag{ mu_key_, tag});
}


void
HydraulicConductivityEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> k = S.GetPtr<CompositeVector>(k_key_, tag);
  Teuchos::RCP<const CompositeVector> rho = S.GetPtr<CompositeVector>(rho_key_, tag);
  Teuchos::RCP<const CompositeVector> mu = S.GetPtr<CompositeVector>(mu_key_, tag);

  const AmanziGeometry::Point& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];

  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    const Epetra_MultiVector& k_v = *k->ViewComponent(*comp, false);
    const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
    const Epetra_MultiVector& mu_v = *mu->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->HydraulicConductivity(k_v[0][i], rho_v[0][i], mu_v[0][i],gz);
    }
  }
}


void
HydraulicConductivityEvaluator::EvaluatePartialDerivative_(const State& S,
        const Key& wrt_key, const Tag& wrt_tag, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> k = S.GetPtr<CompositeVector>(k_key_, tag);
  Teuchos::RCP<const CompositeVector> rho = S.GetPtr<CompositeVector>(rho_key_, tag);
  Teuchos::RCP<const CompositeVector> mu = S.GetPtr<CompositeVector>(mu_key_, tag);

  const AmanziGeometry::Point& gravity = S.Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
  double gz = -gravity[2];

  if (wrt_key == k_key_) {
    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& k_v = *k->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& mu_v = *mu->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DHydraulicConductivityDPermeability(k_v[0][i], rho_v[0][i], mu_v[0][i],gz);
      }
    }

  } else if (wrt_key == rho_key_) {
    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& k_v = *k->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& mu_v = *mu->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DHydraulicConductivityDMassDensityLiquid(k_v[0][i], rho_v[0][i], mu_v[0][i],gz);
      }
    }

  } else if (wrt_key == mu_key_) {
    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& k_v = *k->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& mu_v = *mu->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DHydraulicConductivityDViscosityLiquid(k_v[0][i], rho_v[0][i], mu_v[0][i],gz);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
