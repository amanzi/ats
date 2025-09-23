/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  The surface ice energy evaluator is an algebraic evaluator of a given model.
Energy evaulator for ice+liquid surface water.
  Generated via evaluator_generator.
*/

#include "surface_ice_energy_evaluator.hh"
#include "surface_ice_energy_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {
namespace Relations {

// Constructor from ParameterList
SurfaceIceEnergyEvaluator::SurfaceIceEnergyEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("surface_ice_energy parameters");
  model_ = Teuchos::rcp(new SurfaceIceEnergyModel(sublist));
  InitializeFromPlist_();
}

// Virtual copy constructor
Teuchos::RCP<Evaluator>
SurfaceIceEnergyEvaluator::Clone() const
{
  return Teuchos::rcp(new SurfaceIceEnergyEvaluator(*this));
}


// Initialize by setting up dependencies
void
SurfaceIceEnergyEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: ponded_depth
  h_key_ = Keys::readKey(plist_, domain_name, "ponded depth", "ponded_depth");
  dependencies_.insert(KeyTag{ h_key_, tag });

  // dependency: unfrozen_fraction
  eta_key_ = Keys::readKey(plist_, domain_name, "unfrozen fraction", "unfrozen_fraction");
  dependencies_.insert(KeyTag{ eta_key_, tag });

  // dependency: molar_density_liquid
  nl_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ nl_key_, tag });

  // dependency: internal_energy_liquid
  ul_key_ = Keys::readKey(plist_, domain_name, "internal energy liquid", "internal_energy_liquid");
  dependencies_.insert(KeyTag{ ul_key_, tag });

  // dependency: molar_density_ice
  ni_key_ = Keys::readKey(plist_, domain_name, "molar density ice", "molar_density_ice");
  dependencies_.insert(KeyTag{ ni_key_, tag });

  // dependency: internal_energy_ice
  ui_key_ = Keys::readKey(plist_, domain_name, "internal energy ice", "internal_energy_ice");
  dependencies_.insert(KeyTag{ ui_key_, tag });

  // dependency: cell_volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
}


void
SurfaceIceEnergyEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> h = S.GetPtr<CompositeVector>(h_key_, tag);
  Teuchos::RCP<const CompositeVector> eta = S.GetPtr<CompositeVector>(eta_key_, tag);
  Teuchos::RCP<const CompositeVector> nl = S.GetPtr<CompositeVector>(nl_key_, tag);
  Teuchos::RCP<const CompositeVector> ul = S.GetPtr<CompositeVector>(ul_key_, tag);
  Teuchos::RCP<const CompositeVector> ni = S.GetPtr<CompositeVector>(ni_key_, tag);
  Teuchos::RCP<const CompositeVector> ui = S.GetPtr<CompositeVector>(ui_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
    const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp, false);
    const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
    const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
    const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
    const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = model_->Energy(
        h_v[0][i], eta_v[0][i], nl_v[0][i], ul_v[0][i], ni_v[0][i], ui_v[0][i], cv_v[0][i]);
    }
  }
}


void
SurfaceIceEnergyEvaluator::EvaluatePartialDerivative_(const State& S,
                                                      const Key& wrt_key,
                                                      const Tag& wrt_tag,
                                                      const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> h = S.GetPtr<CompositeVector>(h_key_, tag);
  Teuchos::RCP<const CompositeVector> eta = S.GetPtr<CompositeVector>(eta_key_, tag);
  Teuchos::RCP<const CompositeVector> nl = S.GetPtr<CompositeVector>(nl_key_, tag);
  Teuchos::RCP<const CompositeVector> ul = S.GetPtr<CompositeVector>(ul_key_, tag);
  Teuchos::RCP<const CompositeVector> ni = S.GetPtr<CompositeVector>(ni_key_, tag);
  Teuchos::RCP<const CompositeVector> ui = S.GetPtr<CompositeVector>(ui_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);

  if (wrt_key == h_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDPondedDepth(
          h_v[0][i], eta_v[0][i], nl_v[0][i], ul_v[0][i], ni_v[0][i], ui_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == eta_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDUnfrozenFraction(
          h_v[0][i], eta_v[0][i], nl_v[0][i], ul_v[0][i], ni_v[0][i], ui_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == nl_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDMolarDensityLiquid(
          h_v[0][i], eta_v[0][i], nl_v[0][i], ul_v[0][i], ni_v[0][i], ui_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == ul_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDInternalEnergyLiquid(
          h_v[0][i], eta_v[0][i], nl_v[0][i], ul_v[0][i], ni_v[0][i], ui_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == ni_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDMolarDensityIce(
          h_v[0][i], eta_v[0][i], nl_v[0][i], ul_v[0][i], ni_v[0][i], ui_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == ui_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDInternalEnergyIce(
          h_v[0][i], eta_v[0][i], nl_v[0][i], ul_v[0][i], ni_v[0][i], ui_v[0][i], cv_v[0][i]);
      }
    }

  } else if (wrt_key == cv_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& h_v = *h->ViewComponent(*comp, false);
      const Epetra_MultiVector& eta_v = *eta->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDCellVolume(
          h_v[0][i], eta_v[0][i], nl_v[0][i], ul_v[0][i], ni_v[0][i], ui_v[0][i], cv_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} // namespace Relations
} // namespace Energy
} // namespace ATS_Physics
} // namespace Amanzi
