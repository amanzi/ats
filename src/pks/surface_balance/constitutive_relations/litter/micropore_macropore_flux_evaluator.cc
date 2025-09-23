/*
  Copyright 2010-202x held jointly by participating institutions.
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

.. _evaluator-micropore-macropore-flux-spec:
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
namespace ATS_Physics {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
MicroporeMacroporeFluxEvaluator::MicroporeMacroporeFluxEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("micropore_macropore_flux parameters");
  model_ = Teuchos::rcp(new MicroporeMacroporeFluxModel(sublist));
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
MicroporeMacroporeFluxEvaluator::Clone() const
{
  return Teuchos::rcp(new MicroporeMacroporeFluxEvaluator(*this));
}


// Initialize by setting up dependencies
void
MicroporeMacroporeFluxEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  Key micro_domain = Keys::getDomain(my_keys_.front().first);
  Key macro_domain = Keys::readDomainHint(plist_, micro_domain, "micropore", "macropore");
  auto tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: micropore_pressure
  pm_key_ = Keys::readKey(plist_, micro_domain, "micropore pressure", "pressure");
  dependencies_.insert(KeyTag{ pm_key_, tag });
  // dependency: pressure
  pM_key_ = Keys::readKey(plist_, macro_domain, "macropore pressure", "pressure");
  dependencies_.insert(KeyTag{ pM_key_, tag });
  // dependency: micropore_relative_permeability
  krm_key_ =
    Keys::readKey(plist_, micro_domain, "micropore relative permeability", "relative_permeability");
  dependencies_.insert(KeyTag{ krm_key_, tag });
  // dependency: relative_permeability
  krM_key_ =
    Keys::readKey(plist_, macro_domain, "macropore relative permeability", "relative_permeability");
  dependencies_.insert(KeyTag{ krM_key_, tag });
  // dependency: micropore_absolute_permeability
  K_key_ = Keys::readKey(plist_, macro_domain, "macropore absolute permeability", "permeability");
  dependencies_.insert(KeyTag{ K_key_, tag });
  // dependency: micropore_molar_density_liquid
  den_key_ =
    Keys::readKey(plist_, micro_domain, "micropore molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ den_key_, tag });
}


void
MicroporeMacroporeFluxEvaluator::Evaluate_(const State& S,
                                           const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pm = S.GetPtr<CompositeVector>(pm_key_, tag);
  Teuchos::RCP<const CompositeVector> pM = S.GetPtr<CompositeVector>(pM_key_, tag);
  Teuchos::RCP<const CompositeVector> krM = S.GetPtr<CompositeVector>(krM_key_, tag);
  Teuchos::RCP<const CompositeVector> krm = S.GetPtr<CompositeVector>(krm_key_, tag);
  Teuchos::RCP<const CompositeVector> K = S.GetPtr<CompositeVector>(K_key_, tag);
  Teuchos::RCP<const CompositeVector> den = S.GetPtr<CompositeVector>(den_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
    const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
    const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
    const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
    const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
    const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = den_v[0][i] * model_->MicroporeMacroporeFlux(
                                       pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
    }
  }
}


void
MicroporeMacroporeFluxEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> pm = S.GetPtr<CompositeVector>(pm_key_, tag);
  Teuchos::RCP<const CompositeVector> pM = S.GetPtr<CompositeVector>(pM_key_, tag);
  Teuchos::RCP<const CompositeVector> krM = S.GetPtr<CompositeVector>(krM_key_, tag);
  Teuchos::RCP<const CompositeVector> krm = S.GetPtr<CompositeVector>(krm_key_, tag);
  Teuchos::RCP<const CompositeVector> K = S.GetPtr<CompositeVector>(K_key_, tag);
  Teuchos::RCP<const CompositeVector> den = S.GetPtr<CompositeVector>(den_key_, tag);

  if (wrt_key == pm_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] =
          den_v[0][i] * model_->DMicroporeMacroporeFluxDMicroporePressure(
                          pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == pM_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] =
          den_v[0][i] * model_->DMicroporeMacroporeFluxDPressure(
                          pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == krM_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] =
          den_v[0][i] * model_->DMicroporeMacroporeFluxDRelativePermeability(
                          pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == krm_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] =
          den_v[0][i] * model_->DMicroporeMacroporeFluxDMicroporeRelativePermeability(
                          pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == K_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      const Epetra_MultiVector& den_v = *den->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] =
          den_v[0][i] * model_->DMicroporeMacroporeFluxDMicroporeAbsolutePermeability(
                          pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }
  } else if (wrt_key == den_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& pm_v = *pm->ViewComponent(*comp, false);
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& krm_v = *krm->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->MicroporeMacroporeFlux(
          pm_v[0][i], pM_v[0][i], krM_v[0][i], krm_v[0][i], K_v[0][i]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace Relations
} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
