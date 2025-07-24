/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

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
RichardsEnergyEvaluator::RichardsEnergyEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("richards_energy parameters");
  model_ = Teuchos::rcp(new RichardsEnergyModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
RichardsEnergyEvaluator::RichardsEnergyEvaluator(const RichardsEnergyEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other),
    phi_key_(other.phi_key_),
    phi0_key_(other.phi0_key_),
    sl_key_(other.sl_key_),
    nl_key_(other.nl_key_),
    ul_key_(other.ul_key_),
    rho_r_key_(other.rho_r_key_),
    ur_key_(other.ur_key_),
    cv_key_(other.cv_key_),
    model_(other.model_)
{}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
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
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: porosity
  phi_key_ = Keys::readKey(plist_, domain_name, "porosity", "porosity");
  dependencies_.insert(KeyTag{ phi_key_, tag });

  // dependency: base_porosity
  phi0_key_ = Keys::readKey(plist_, domain_name, "base porosity", "base_porosity");
  dependencies_.insert(KeyTag{ phi0_key_, tag });

  // dependency: saturation_liquid
  sl_key_ = Keys::readKey(plist_, domain_name, "saturation liquid", "saturation_liquid");
  dependencies_.insert(KeyTag{ sl_key_, tag });

  // dependency: molar_density_liquid
  nl_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ nl_key_, tag });

  // dependency: internal_energy_liquid
  ul_key_ = Keys::readKey(plist_, domain_name, "internal energy liquid", "internal_energy_liquid");
  dependencies_.insert(KeyTag{ ul_key_, tag });

  // dependency: density_rock
  rho_r_key_ = Keys::readKey(plist_, domain_name, "density rock", "density_rock");
  dependencies_.insert(KeyTag{ rho_r_key_, tag });

  // dependency: internal_energy_rock
  ur_key_ = Keys::readKey(plist_, domain_name, "internal energy rock", "internal_energy_rock");
  dependencies_.insert(KeyTag{ ur_key_, tag });

  // dependency: cell_volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
}


void
RichardsEnergyEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(phi_key_, tag);
  Teuchos::RCP<const CompositeVector> phi0 = S.GetPtr<CompositeVector>(phi0_key_, tag);
  Teuchos::RCP<const CompositeVector> sl = S.GetPtr<CompositeVector>(sl_key_, tag);
  Teuchos::RCP<const CompositeVector> nl = S.GetPtr<CompositeVector>(nl_key_, tag);
  Teuchos::RCP<const CompositeVector> ul = S.GetPtr<CompositeVector>(ul_key_, tag);
  Teuchos::RCP<const CompositeVector> rho_r = S.GetPtr<CompositeVector>(rho_r_key_, tag);
  Teuchos::RCP<const CompositeVector> ur = S.GetPtr<CompositeVector>(ur_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
    const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
    const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
    const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
    const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
    const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = model_->Energy(phi_v[0][i],
                                      phi0_v[0][i],
                                      sl_v[0][i],
                                      nl_v[0][i],
                                      ul_v[0][i],
                                      rho_r_v[0][i],
                                      ur_v[0][i],
                                      cv_v[0][i]);
    }
  }
}


void
RichardsEnergyEvaluator::EvaluatePartialDerivative_(const State& S,
                                                    const Key& wrt_key,
                                                    const Tag& wrt_tag,
                                                    const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(phi_key_, tag);
  Teuchos::RCP<const CompositeVector> phi0 = S.GetPtr<CompositeVector>(phi0_key_, tag);
  Teuchos::RCP<const CompositeVector> sl = S.GetPtr<CompositeVector>(sl_key_, tag);
  Teuchos::RCP<const CompositeVector> nl = S.GetPtr<CompositeVector>(nl_key_, tag);
  Teuchos::RCP<const CompositeVector> ul = S.GetPtr<CompositeVector>(ul_key_, tag);
  Teuchos::RCP<const CompositeVector> rho_r = S.GetPtr<CompositeVector>(rho_r_key_, tag);
  Teuchos::RCP<const CompositeVector> ur = S.GetPtr<CompositeVector>(ur_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);

  if (wrt_key == phi_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDPorosity(phi_v[0][i],
                                                  phi0_v[0][i],
                                                  sl_v[0][i],
                                                  nl_v[0][i],
                                                  ul_v[0][i],
                                                  rho_r_v[0][i],
                                                  ur_v[0][i],
                                                  cv_v[0][i]);
      }
    }

  } else if (wrt_key == phi0_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDBasePorosity(phi_v[0][i],
                                                      phi0_v[0][i],
                                                      sl_v[0][i],
                                                      nl_v[0][i],
                                                      ul_v[0][i],
                                                      rho_r_v[0][i],
                                                      ur_v[0][i],
                                                      cv_v[0][i]);
      }
    }

  } else if (wrt_key == sl_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDSaturationLiquid(phi_v[0][i],
                                                          phi0_v[0][i],
                                                          sl_v[0][i],
                                                          nl_v[0][i],
                                                          ul_v[0][i],
                                                          rho_r_v[0][i],
                                                          ur_v[0][i],
                                                          cv_v[0][i]);
      }
    }

  } else if (wrt_key == nl_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDMolarDensityLiquid(phi_v[0][i],
                                                            phi0_v[0][i],
                                                            sl_v[0][i],
                                                            nl_v[0][i],
                                                            ul_v[0][i],
                                                            rho_r_v[0][i],
                                                            ur_v[0][i],
                                                            cv_v[0][i]);
      }
    }

  } else if (wrt_key == ul_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDInternalEnergyLiquid(phi_v[0][i],
                                                              phi0_v[0][i],
                                                              sl_v[0][i],
                                                              nl_v[0][i],
                                                              ul_v[0][i],
                                                              rho_r_v[0][i],
                                                              ur_v[0][i],
                                                              cv_v[0][i]);
      }
    }

  } else if (wrt_key == rho_r_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDDensityRock(phi_v[0][i],
                                                     phi0_v[0][i],
                                                     sl_v[0][i],
                                                     nl_v[0][i],
                                                     ul_v[0][i],
                                                     rho_r_v[0][i],
                                                     ur_v[0][i],
                                                     cv_v[0][i]);
      }
    }

  } else if (wrt_key == ur_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDInternalEnergyRock(phi_v[0][i],
                                                            phi0_v[0][i],
                                                            sl_v[0][i],
                                                            nl_v[0][i],
                                                            ul_v[0][i],
                                                            rho_r_v[0][i],
                                                            ur_v[0][i],
                                                            cv_v[0][i]);
      }
    }

  } else if (wrt_key == cv_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDCellVolume(phi_v[0][i],
                                                    phi0_v[0][i],
                                                    sl_v[0][i],
                                                    nl_v[0][i],
                                                    ul_v[0][i],
                                                    rho_r_v[0][i],
                                                    ur_v[0][i],
                                                    cv_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace Relations
} // namespace Energy
} // namespace Amanzi
