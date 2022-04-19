/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! Exchange flux between multiple continua.

/*!

Evaluates the following exchange flux model:

.. math::
   q_{exchange} = n_l k_r K \frac{\Gamma}{\delta} (p_M - p_m)

where :math:`p` is the pressure of the Macro and micro porespaces,
respectively, K is some measure of an absolute permeability, :math:`\Gamma [-]`
is the exchange coefficient, :math:`\delta [m]` is a unit of distance
characterizing the typical distance between pore, and :math:`k_r` is the
relative permeability, which is upwinded based on the larger of the two
pressures, 'n_l' is molar density liquid 

Note that the expected domain for this is the micropore domain, but may be
changed on the input line.

.. _micropore-macropore-flux-evaluator-spec:
.. admonition:: micropore-micropore-flux-evaluator

   * `"micropore domain`" ``[string]`` **""** Defaults to the domain of the flux's
     variable name.

   * `"macropore domain`" ``[string]`` **macropore**

   * `"micropore macropore flux model parameters`" ``[micropore-macropore-flux-model-spec]``

   KEYS:

   * `"micropore pressure`" **pressure**
   * `"macropore pressure`" **MACROPORE_DOMAIN-pressure**
   * `"micropore relative permeability`" **relative_permeability**
   * `"macropore relative permeability`" **MACROPORE_DOMAIN-relative_permeability**
   * `"permeability`" **permeability**
   * `"micropore molar density liquid`" **molar_density_liquid** 

*/

#include "micropore_macropore_flux_evaluator.hh"
#include "micropore_macropore_flux_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
MicroporeMacroporeFluxEvaluator::MicroporeMacroporeFluxEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("micropore_macropore_flux parameters");
  model_ = Teuchos::rcp(new MicroporeMacroporeFluxModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
MicroporeMacroporeFluxEvaluator::MicroporeMacroporeFluxEvaluator(const MicroporeMacroporeFluxEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pm_key_(other.pm_key_),
    pM_key_(other.pM_key_),
    krM_key_(other.krM_key_),
    krm_key_(other.krm_key_),
    K_key_(other.K_key_),
    den_key_(other.den_key_),
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
MicroporeMacroporeFluxEvaluator::Clone() const
{
  return Teuchos::rcp(new MicroporeMacroporeFluxEvaluator(*this));
}


// Initialize by setting up dependencies
void
MicroporeMacroporeFluxEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  //Key micro_domain = Keys::getDomain(my_key_);
  Key micro_domain = plist_.get<std::string>("micropore domain name", "domain");
  Key macro_domain = plist_.get<std::string>("macropore domain name", "macropore");

  // - pull Keys from plist
  // dependency: micropore_pressure
  pm_key_ = Keys::readKey(plist_, micro_domain, "micropore pressure", "pressure");
  dependencies_.insert(pm_key_);
  // dependency: pressure
  pM_key_ = Keys::readKey(plist_, macro_domain, "macropore pressure", "pressure");
  dependencies_.insert(pM_key_);
  // dependency: micropore_relative_permeability
  krm_key_ = Keys::readKey(plist_, micro_domain, "micropore relative permeability", "relative_permeability");
  dependencies_.insert(krm_key_);
  // dependency: relative_permeability
  krM_key_ = Keys::readKey(plist_, macro_domain, "macropore relative permeability", "relative_permeability");
  dependencies_.insert(krM_key_);
  // dependency: micropore_absolute_permeability
  K_key_ = Keys::readKey(plist_, macro_domain, "macropore absolute permeability", "permeability");  
  dependencies_.insert(K_key_);
  // dependency: micropore_molar_density_liquid
  den_key_ = Keys::readKey(plist_, micro_domain, "micropore molar density liquid", "molar_density_liquid");
  dependencies_.insert(den_key_);
    
}


void
MicroporeMacroporeFluxEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> pm = S->GetFieldData(pm_key_);
  Teuchos::RCP<const CompositeVector> pM = S->GetFieldData(pM_key_);
  Teuchos::RCP<const CompositeVector> krM = S->GetFieldData(krM_key_);
  Teuchos::RCP<const CompositeVector> krm = S->GetFieldData(krm_key_);
  Teuchos::RCP<const CompositeVector> K = S->GetFieldData(K_key_);
  Teuchos::RCP<const CompositeVector> den = S->GetFieldData(den_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
    const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
    const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
    const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
    const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
    const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = den_v[0][i]*model_->MicroporeMacroporeFlux(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
    }
  }
}


void
MicroporeMacroporeFluxEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> pm = S->GetFieldData(pm_key_);
  Teuchos::RCP<const CompositeVector> pM = S->GetFieldData(pM_key_);
  Teuchos::RCP<const CompositeVector> krM = S->GetFieldData(krM_key_);
  Teuchos::RCP<const CompositeVector> krm = S->GetFieldData(krm_key_);
  Teuchos::RCP<const CompositeVector> K = S->GetFieldData(K_key_);
  Teuchos::RCP<const CompositeVector> den = S->GetFieldData(den_key_);

  if (wrt_key == pm_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = den_v[0][i]*model_->DMicroporeMacroporeFluxDMicroporePressure(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == pM_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = den_v[0][i]*model_->DMicroporeMacroporeFluxDPressure(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == krM_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = den_v[0][i]*model_->DMicroporeMacroporeFluxDRelativePermeability(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == krm_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = den_v[0][i]*model_->DMicroporeMacroporeFluxDMicroporeRelativePermeability(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == K_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = den_v[0][i]*model_->DMicroporeMacroporeFluxDMicroporeAbsolutePermeability(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }
  } else if (wrt_key == den_key_){
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->MicroporeMacroporeFlux(pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

} //namespace
} //namespace
} //namespace
