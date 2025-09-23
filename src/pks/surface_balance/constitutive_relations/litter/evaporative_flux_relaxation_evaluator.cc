/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The evaporative flux relaxation evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    myKeyFirst = evaporative
    evalNameString = evaporative flux relaxation
    evalNameCaps = EVAPORATIVE_FLUX_RELAXATION
    namespaceCaps = SURFACEBALANCE
    evalClassName = EvaporativeFluxRelaxation
    namespace = SurfaceBalance
    myMethodDeclarationArgs = double wc, double rho, double L
    myKey = evaporative_flux
    evalName = evaporative_flux_relaxation
    modelMethodDeclaration =   double EvaporativeFlux(double wc, double rho, double L) const;
    myKeyMethod = EvaporativeFlux
    myMethodArgs = wc_v[0][i], rho_v[0][i], L_v[0][i]

*/

#include "evaporative_flux_relaxation_evaluator.hh"
#include "evaporative_flux_relaxation_model.hh"

namespace Amanzi {
namespace ATS_Physics {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
EvaporativeFluxRelaxationEvaluator::EvaporativeFluxRelaxationEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("evaporative_flux_relaxation parameters");
  model_ = Teuchos::rcp(new EvaporativeFluxRelaxationModel(sublist));
  InitializeFromPlist_();
}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
EvaporativeFluxRelaxationEvaluator::Clone() const
{
  return Teuchos::rcp(new EvaporativeFluxRelaxationEvaluator(*this));
}


// Initialize by setting up dependencies
void
EvaporativeFluxRelaxationEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  std::string domain_name = Keys::getDomain(my_keys_.front().first);
  auto tag = my_keys_.front().second;

  // - pull Keys from plist
  // dependency: litter_water_content
  wc_key_ = Keys::readKey(plist_, domain_name, "litter water content", "water_content");
  dependencies_.insert(KeyTag{ wc_key_, tag });

  // dependency: surface_molar_density_liquid
  rho_key_ = Keys::readKey(plist_, domain_name, "molar density liquid", "molar_density_liquid");
  dependencies_.insert(KeyTag{ rho_key_, tag });

  // dependency: litter_thickness
  thickness_key_ = Keys::readKey(plist_, domain_name, "litter thickness", "thickness");
  dependencies_.insert(KeyTag{ thickness_key_, tag });

  // dependency: cell volume
  cv_key_ = Keys::readKey(plist_, domain_name, "cell volume", "cell_volume");
  dependencies_.insert(KeyTag{ cv_key_, tag });
}


void
EvaporativeFluxRelaxationEvaluator::Evaluate_(const State& S,
                                              const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> wc = S.GetPtr<CompositeVector>(wc_key_, tag);
  Teuchos::RCP<const CompositeVector> rho = S.GetPtr<CompositeVector>(rho_key_, tag);
  Teuchos::RCP<const CompositeVector> L = S.GetPtr<CompositeVector>(thickness_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp, false);
    const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
    const Epetra_MultiVector& L_v = *L->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    int ncomp = result[0]->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = model_->EvaporativeFlux(wc_v[0][i] / cv_v[0][i], rho_v[0][i], L_v[0][i]);
    }
  }
}


void
EvaporativeFluxRelaxationEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  auto tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> wc = S.GetPtr<CompositeVector>(wc_key_, tag);
  Teuchos::RCP<const CompositeVector> rho = S.GetPtr<CompositeVector>(rho_key_, tag);
  Teuchos::RCP<const CompositeVector> L = S.GetPtr<CompositeVector>(thickness_key_, tag);
  Teuchos::RCP<const CompositeVector> cv = S.GetPtr<CompositeVector>(cv_key_, tag);

  if (wrt_key == wc_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& L_v = *L->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEvaporativeFluxDLitterWaterContent(
                           wc_v[0][i] / cv_v[0][i], rho_v[0][i], L_v[0][i]) /
                         cv_v[0][i];
      }
    }

  } else if (wrt_key == rho_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& L_v = *L->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DEvaporativeFluxDSurfaceMolarDensityLiquid(
          wc_v[0][i] / cv_v[0][i], rho_v[0][i], L_v[0][i]);
      }
    }

  } else if (wrt_key == thickness_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin() ; comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& L_v = *L->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

      int ncomp = result[0]->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] =
          model_->DEvaporativeFluxDLitterThickness(wc_v[0][i] / cv_v[0][i], rho_v[0][i], L_v[0][i]);
      }
    }
  }
}


} // namespace Relations
} // namespace SurfaceBalance
} // namespace ATS_Physics
} // namespace Amanzi
