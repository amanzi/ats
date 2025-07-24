/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  The liquid + gas water content evaluator is an algebraic evaluator of a given model.
Water content for a two-phase, liquid+water vapor evaluator.
  Generated via evaluator_generator.
*/

#include "liquid_gas_water_content_evaluator.hh"
#include "liquid_gas_water_content_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

// Constructor from ParameterList
LiquidGasWaterContentEvaluator::LiquidGasWaterContentEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("liquid_gas_water_content parameters");
  model_ = Teuchos::rcp(new LiquidGasWaterContentModel(sublist));
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
LiquidGasWaterContentEvaluator::Clone() const
{
  return Teuchos::rcp(new LiquidGasWaterContentEvaluator(*this));
}


// Initialize by setting up dependencies
void
LiquidGasWaterContentEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: porosity
  phi_key_ = Keys::readKey(plist_, domain_name, "porosity", "porosity");
  dependencies_.insert(KeyTag{ phi_key_, tag });

  // dependency: saturation_liquid
  sl_key_ = Keys::readKey(plist_, domain_name, "saturation liquid", "saturation_liquid");
  dependencies_.insert(KeyTag{ sl_key_, tag });

  // dependency: molar_density_liquid
  nl_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ nl_key_, tag });

  // dependency: saturation_gas
  sg_key_ = Keys::readKey(plist_, domain_name, "saturation gas", "saturation_gas");
  dependencies_.insert(KeyTag{ sg_key_, tag });

  // dependency: molar_density_gas
  ng_key_ = Keys::readKey(plist_, domain_name, "molar density gas", "molar_density_gas");
  dependencies_.insert(KeyTag{ ng_key_, tag });

  // dependency: mol_frac_gas
  omega_key_ = Keys::readKey(plist_, domain_name, "mol frac gas", "mol_frac_gas");
  dependencies_.insert(KeyTag{ omega_key_, tag });

  // dependency: cell_volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
}


void
LiquidGasWaterContentEvaluator::Evaluate_(const State& S,
                                          const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(phi_key_, tag);
  Teuchos::RCP<const CompositeVector> sl = S.GetPtr<CompositeVector>(sl_key_, tag);
  Teuchos::RCP<const CompositeVector> nl = S.GetPtr<CompositeVector>(nl_key_, tag);
  Teuchos::RCP<const CompositeVector> sg = S.GetPtr<CompositeVector>(sg_key_, tag);
  Teuchos::RCP<const CompositeVector> ng = S.GetPtr<CompositeVector>(ng_key_, tag);
  Teuchos::RCP<const CompositeVector> omega = S.GetPtr<CompositeVector>(omega_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
    const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
    const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
    const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
    const Epetra_MultiVector& omega_v = *omega->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = model_->WaterContent(
        phi_v[0][i], sl_v[0][i], nl_v[0][i], sg_v[0][i], ng_v[0][i], omega_v[0][i], cv_v[0][i]);
    }
  }
}


void
LiquidGasWaterContentEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(phi_key_, tag);
  Teuchos::RCP<const CompositeVector> sl = S.GetPtr<CompositeVector>(sl_key_, tag);
  Teuchos::RCP<const CompositeVector> nl = S.GetPtr<CompositeVector>(nl_key_, tag);
  Teuchos::RCP<const CompositeVector> sg = S.GetPtr<CompositeVector>(sg_key_, tag);
  Teuchos::RCP<const CompositeVector> ng = S.GetPtr<CompositeVector>(ng_key_, tag);
  Teuchos::RCP<const CompositeVector> omega = S.GetPtr<CompositeVector>(omega_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);

  if (wrt_key == phi_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& omega_v = *omega->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDPorosity(
          phi_v[0][i], sl_v[0][i], nl_v[0][i], sg_v[0][i], ng_v[0][i], omega_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == sl_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& omega_v = *omega->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDSaturationLiquid(
          phi_v[0][i], sl_v[0][i], nl_v[0][i], sg_v[0][i], ng_v[0][i], omega_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == nl_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& omega_v = *omega->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDMolarDensityLiquid(
          phi_v[0][i], sl_v[0][i], nl_v[0][i], sg_v[0][i], ng_v[0][i], omega_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == sg_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& omega_v = *omega->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDSaturationGas(
          phi_v[0][i], sl_v[0][i], nl_v[0][i], sg_v[0][i], ng_v[0][i], omega_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == ng_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& omega_v = *omega->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDMolarDensityGas(
          phi_v[0][i], sl_v[0][i], nl_v[0][i], sg_v[0][i], ng_v[0][i], omega_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == omega_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& omega_v = *omega->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDMolFracGas(
          phi_v[0][i], sl_v[0][i], nl_v[0][i], sg_v[0][i], ng_v[0][i], omega_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == cv_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& omega_v = *omega->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DWaterContentDCellVolume(
          phi_v[0][i], sl_v[0][i], nl_v[0][i], sg_v[0][i], ng_v[0][i], omega_v[0][i], cv_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
