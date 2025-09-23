/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  The three phase energy evaluator is an algebraic evaluator of a given model.
Energy for a three-phase, gas+liquid+ice evaluator.
  Generated via evaluator_generator.
*/

#include "three_phase_energy_evaluator.hh"
#include "three_phase_energy_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace Energy {
namespace Relations {

// Constructor from ParameterList
ThreePhaseEnergyEvaluator::ThreePhaseEnergyEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("three_phase_energy parameters");
  model_ = Teuchos::rcp(new ThreePhaseEnergyModel(sublist));
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
ThreePhaseEnergyEvaluator::Clone() const
{
  return Teuchos::rcp(new ThreePhaseEnergyEvaluator(*this));
}


// Initialize by setting up dependencies
void
ThreePhaseEnergyEvaluator::InitializeFromPlist_()
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

  // dependency: saturation_ice
  si_key_ = Keys::readKey(plist_, domain_name, "saturation ice", "saturation_ice");
  dependencies_.insert(KeyTag{ si_key_, tag });

  // dependency: molar_density_ice
  ni_key_ = Keys::readKey(plist_, domain_name, "molar density ice", "molar_density_ice");
  dependencies_.insert(KeyTag{ ni_key_, tag });

  // dependency: internal_energy_ice
  ui_key_ = Keys::readKey(plist_, domain_name, "internal energy ice", "internal_energy_ice");
  dependencies_.insert(KeyTag{ ui_key_, tag });

  // dependency: saturation_gas
  sg_key_ = Keys::readKey(plist_, domain_name, "saturation gas", "saturation_gas");
  dependencies_.insert(KeyTag{ sg_key_, tag });

  // dependency: molar_density_gas
  ng_key_ = Keys::readKey(plist_, domain_name, "molar density gas", "molar_density_gas");
  dependencies_.insert(KeyTag{ ng_key_, tag });

  // dependency: internal_energy_gas
  ug_key_ = Keys::readKey(plist_, domain_name, "internal energy gas", "internal_energy_gas");
  dependencies_.insert(KeyTag{ ug_key_, tag });

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
ThreePhaseEnergyEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> phi = S.GetPtr<CompositeVector>(phi_key_, tag);
  Teuchos::RCP<const CompositeVector> phi0 = S.GetPtr<CompositeVector>(phi0_key_, tag);
  Teuchos::RCP<const CompositeVector> sl = S.GetPtr<CompositeVector>(sl_key_, tag);
  Teuchos::RCP<const CompositeVector> nl = S.GetPtr<CompositeVector>(nl_key_, tag);
  Teuchos::RCP<const CompositeVector> ul = S.GetPtr<CompositeVector>(ul_key_, tag);
  Teuchos::RCP<const CompositeVector> si = S.GetPtr<CompositeVector>(si_key_, tag);
  Teuchos::RCP<const CompositeVector> ni = S.GetPtr<CompositeVector>(ni_key_, tag);
  Teuchos::RCP<const CompositeVector> ui = S.GetPtr<CompositeVector>(ui_key_, tag);
  Teuchos::RCP<const CompositeVector> sg = S.GetPtr<CompositeVector>(sg_key_, tag);
  Teuchos::RCP<const CompositeVector> ng = S.GetPtr<CompositeVector>(ng_key_, tag);
  Teuchos::RCP<const CompositeVector> ug = S.GetPtr<CompositeVector>(ug_key_, tag);
  Teuchos::RCP<const CompositeVector> rho_r = S.GetPtr<CompositeVector>(rho_r_key_, tag);
  Teuchos::RCP<const CompositeVector> ur = S.GetPtr<CompositeVector>(ur_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
    const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
    const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
    const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
    const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
    const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
    const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
    const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
    const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
    const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
    const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
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
                                      si_v[0][i],
                                      ni_v[0][i],
                                      ui_v[0][i],
                                      sg_v[0][i],
                                      ng_v[0][i],
                                      ug_v[0][i],
                                      rho_r_v[0][i],
                                      ur_v[0][i],
                                      cv_v[0][i]);
    }
  }
}


void
ThreePhaseEnergyEvaluator::EvaluatePartialDerivative_(const State& S,
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
  Teuchos::RCP<const CompositeVector> si = S.GetPtr<CompositeVector>(si_key_, tag);
  Teuchos::RCP<const CompositeVector> ni = S.GetPtr<CompositeVector>(ni_key_, tag);
  Teuchos::RCP<const CompositeVector> ui = S.GetPtr<CompositeVector>(ui_key_, tag);
  Teuchos::RCP<const CompositeVector> sg = S.GetPtr<CompositeVector>(sg_key_, tag);
  Teuchos::RCP<const CompositeVector> ng = S.GetPtr<CompositeVector>(ng_key_, tag);
  Teuchos::RCP<const CompositeVector> ug = S.GetPtr<CompositeVector>(ug_key_, tag);
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
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
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
                                                  si_v[0][i],
                                                  ni_v[0][i],
                                                  ui_v[0][i],
                                                  sg_v[0][i],
                                                  ng_v[0][i],
                                                  ug_v[0][i],
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
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
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
                                                      si_v[0][i],
                                                      ni_v[0][i],
                                                      ui_v[0][i],
                                                      sg_v[0][i],
                                                      ng_v[0][i],
                                                      ug_v[0][i],
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
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
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
                                                          si_v[0][i],
                                                          ni_v[0][i],
                                                          ui_v[0][i],
                                                          sg_v[0][i],
                                                          ng_v[0][i],
                                                          ug_v[0][i],
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
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
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
                                                            si_v[0][i],
                                                            ni_v[0][i],
                                                            ui_v[0][i],
                                                            sg_v[0][i],
                                                            ng_v[0][i],
                                                            ug_v[0][i],
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
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
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
                                                              si_v[0][i],
                                                              ni_v[0][i],
                                                              ui_v[0][i],
                                                              sg_v[0][i],
                                                              ng_v[0][i],
                                                              ug_v[0][i],
                                                              rho_r_v[0][i],
                                                              ur_v[0][i],
                                                              cv_v[0][i]);
      }
    }

  } else if (wrt_key == si_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDSaturationIce(phi_v[0][i],
                                                       phi0_v[0][i],
                                                       sl_v[0][i],
                                                       nl_v[0][i],
                                                       ul_v[0][i],
                                                       si_v[0][i],
                                                       ni_v[0][i],
                                                       ui_v[0][i],
                                                       sg_v[0][i],
                                                       ng_v[0][i],
                                                       ug_v[0][i],
                                                       rho_r_v[0][i],
                                                       ur_v[0][i],
                                                       cv_v[0][i]);
      }
    }

  } else if (wrt_key == ni_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDMolarDensityIce(phi_v[0][i],
                                                         phi0_v[0][i],
                                                         sl_v[0][i],
                                                         nl_v[0][i],
                                                         ul_v[0][i],
                                                         si_v[0][i],
                                                         ni_v[0][i],
                                                         ui_v[0][i],
                                                         sg_v[0][i],
                                                         ng_v[0][i],
                                                         ug_v[0][i],
                                                         rho_r_v[0][i],
                                                         ur_v[0][i],
                                                         cv_v[0][i]);
      }
    }

  } else if (wrt_key == ui_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDInternalEnergyIce(phi_v[0][i],
                                                           phi0_v[0][i],
                                                           sl_v[0][i],
                                                           nl_v[0][i],
                                                           ul_v[0][i],
                                                           si_v[0][i],
                                                           ni_v[0][i],
                                                           ui_v[0][i],
                                                           sg_v[0][i],
                                                           ng_v[0][i],
                                                           ug_v[0][i],
                                                           rho_r_v[0][i],
                                                           ur_v[0][i],
                                                           cv_v[0][i]);
      }
    }

  } else if (wrt_key == sg_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDSaturationGas(phi_v[0][i],
                                                       phi0_v[0][i],
                                                       sl_v[0][i],
                                                       nl_v[0][i],
                                                       ul_v[0][i],
                                                       si_v[0][i],
                                                       ni_v[0][i],
                                                       ui_v[0][i],
                                                       sg_v[0][i],
                                                       ng_v[0][i],
                                                       ug_v[0][i],
                                                       rho_r_v[0][i],
                                                       ur_v[0][i],
                                                       cv_v[0][i]);
      }
    }

  } else if (wrt_key == ng_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDMolarDensityGas(phi_v[0][i],
                                                         phi0_v[0][i],
                                                         sl_v[0][i],
                                                         nl_v[0][i],
                                                         ul_v[0][i],
                                                         si_v[0][i],
                                                         ni_v[0][i],
                                                         ui_v[0][i],
                                                         sg_v[0][i],
                                                         ng_v[0][i],
                                                         ug_v[0][i],
                                                         rho_r_v[0][i],
                                                         ur_v[0][i],
                                                         cv_v[0][i]);
      }
    }

  } else if (wrt_key == ug_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& phi_v = *phi->ViewComponent(*comp, false);
      const Epetra_MultiVector& phi0_v = *phi0->ViewComponent(*comp, false);
      const Epetra_MultiVector& sl_v = *sl->ViewComponent(*comp, false);
      const Epetra_MultiVector& nl_v = *nl->ViewComponent(*comp, false);
      const Epetra_MultiVector& ul_v = *ul->ViewComponent(*comp, false);
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_r_v = *rho_r->ViewComponent(*comp, false);
      const Epetra_MultiVector& ur_v = *ur->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEnergyDInternalEnergyGas(phi_v[0][i],
                                                           phi0_v[0][i],
                                                           sl_v[0][i],
                                                           nl_v[0][i],
                                                           ul_v[0][i],
                                                           si_v[0][i],
                                                           ni_v[0][i],
                                                           ui_v[0][i],
                                                           sg_v[0][i],
                                                           ng_v[0][i],
                                                           ug_v[0][i],
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
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
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
                                                     si_v[0][i],
                                                     ni_v[0][i],
                                                     ui_v[0][i],
                                                     sg_v[0][i],
                                                     ng_v[0][i],
                                                     ug_v[0][i],
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
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
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
                                                            si_v[0][i],
                                                            ni_v[0][i],
                                                            ui_v[0][i],
                                                            sg_v[0][i],
                                                            ng_v[0][i],
                                                            ug_v[0][i],
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
      const Epetra_MultiVector& si_v = *si->ViewComponent(*comp, false);
      const Epetra_MultiVector& ni_v = *ni->ViewComponent(*comp, false);
      const Epetra_MultiVector& ui_v = *ui->ViewComponent(*comp, false);
      const Epetra_MultiVector& sg_v = *sg->ViewComponent(*comp, false);
      const Epetra_MultiVector& ng_v = *ng->ViewComponent(*comp, false);
      const Epetra_MultiVector& ug_v = *ug->ViewComponent(*comp, false);
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
                                                    si_v[0][i],
                                                    ni_v[0][i],
                                                    ui_v[0][i],
                                                    sg_v[0][i],
                                                    ng_v[0][i],
                                                    ug_v[0][i],
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
} // namespace ATS_Physics
} // namespace Amanzi
