/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/* -----------------------------------------------------------------------------
ATS

Evaluator for water content.

INTERFROST's comparison uses a very odd compressibility term that doesn't
quite fit into either compressible porosity or into a compressible density, so
it needs a special evaluator.
----------------------------------------------------------------------------- */


#include "interfrost_water_content.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

InterfrostWaterContent::InterfrostWaterContent(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // dependency: porosity
  phi_key_ = Keys::readKey(plist_, domain_name, "porosity", "porosity");
  dependencies_.insert(KeyTag{ phi_key_, tag });

  // dependency: saturation_liquid
  sl_key_ = Keys::readKey(plist_, domain_name, "saturation liquid", "saturation_liquid");
  dependencies_.insert(KeyTag{ sl_key_, tag });

  // dependency: molar_density_liquid
  nl_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ nl_key_, tag });

  // dependency: saturation_ice
  si_key_ = Keys::readKey(plist_, domain_name, "saturation ice", "saturation_ice");
  dependencies_.insert(KeyTag{ sl_key_, tag });

  // dependency: molar_density_ice
  ni_key_ = Keys::readKey(plist_, domain_name, "molar density ice", "molar_density_ice");
  dependencies_.insert(KeyTag{ ni_key_, tag });

  // dependency: cell_volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  // dependency: cell_volume
  pres_key_ = Keys::readKey(plist_, domain_name, "pressure", "pressure");
  dependencies_.insert(KeyTag{ cv_key_, tag });

  beta_ = plist.get<double>("compressibility [1/Pa]");
};

Teuchos::RCP<Evaluator>
InterfrostWaterContent::Clone() const
{
  return Teuchos::rcp(new InterfrostWaterContent(*this));
};


void
InterfrostWaterContent::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(phi_key_, tag);
  Teuchos::RCP<const CompositeVector> sl = S.GetPtr<CompositeVector>(sl_key_, tag);
  Teuchos::RCP<const CompositeVector> nl = S.GetPtr<CompositeVector>(nl_key_, tag);
  Teuchos::RCP<const CompositeVector> si = S.GetPtr<CompositeVector>(si_key_, tag);
  Teuchos::RCP<const CompositeVector> ni = S.GetPtr<CompositeVector>(ni_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);
  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
    const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
    const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
    const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      double pc = std::max(pres_v[0][i] - 101325., 0.);
      result_v[0][i] = phi_v[0][i] * (sl_v[0][i] * nl_v[0][i] * (1 + beta_ * pc) + si_v[0][i] * ni_v[0][i]);
      result_v[0][i] *= cv_v[0][i];
    }
  }
};


void
InterfrostWaterContent::EvaluatePartialDerivative_(const State& S,
                                                   const Key& wrt_key,
                                                   const Tag& wrt_tag,
                                                   const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(phi_key_, tag);
  Teuchos::RCP<const CompositeVector> sl = S.GetPtr<CompositeVector>(sl_key_, tag);
  Teuchos::RCP<const CompositeVector> nl = S.GetPtr<CompositeVector>(nl_key_, tag);
  Teuchos::RCP<const CompositeVector> si = S.GetPtr<CompositeVector>(si_key_, tag);
  Teuchos::RCP<const CompositeVector> ni = S.GetPtr<CompositeVector>(ni_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);
  Teuchos::RCP<const CompositeVector> pres = S.GetPtr<CompositeVector>(pres_key_, tag);

  if (wrt_key == phi_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        double pc = std::max(pres_v[0][i] - 101325., 0.);
        result_v[0][i] = (sl_v[0][i] * nl_v[0][i] * (1 + beta_ * pc) + si_v[0][i] * ni_v[0][i]) * cv_v[0][i];
      }
    }
  } else if (wrt_key == sl_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        double pc = std::max(pres_v[0][i] - 101325., 0.);
        result_v[0][i] = nl_v[0][i] * (1 + beta_ * pc) * cv_v[0][i] * phi_v[0][i];
      }
    }
  } else if (wrt_key == nl_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        double pc = std::max(pres_v[0][i] - 101325., 0.);
        result_v[0][i] = sl_v[0][i] * (1 + beta_ * pc) * cv_v[0][i] * phi_v[0][i];
      }
    }
  } else if (wrt_key == si_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = ni_v[0][i] * cv_v[0][i] * phi_v[0][i];
      }
    }
  } else if (wrt_key == ni_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = si_v[0][i] * cv_v[0][i] * phi_v[0][i];
      }
    }
  } else if (wrt_key == pres_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = beta_ * sl_v[0][i] * nl_v[0][i] * cv_v[0][i] * phi_v[0][i];
      }
    }
  } else if (wrt_key == cv_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        double pc = std::max(pres_v[0][i] - 101325., 0.);
        result_v[0][i] = phi_v[0][i] * (sl_v[0][i] * nl_v[0][i] * (1 + beta_ * pc) + si_v[0][i] * ni_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} // namespace Relations
} // namespace Flow
} // namespace Amanzi
