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
  
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "evaporative_flux_relaxation_evaluator.hh"
#include "evaporative_flux_relaxation_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
EvaporativeFluxRelaxationEvaluator::EvaporativeFluxRelaxationEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("evaporative_flux_relaxation parameters");
  model_ = Teuchos::rcp(new EvaporativeFluxRelaxationModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
EvaporativeFluxRelaxationEvaluator::EvaporativeFluxRelaxationEvaluator(const EvaporativeFluxRelaxationEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wc_key_(other.wc_key_),
    rho_key_(other.rho_key_),
    L_key_(other.L_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
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
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);

  std::string my_key_first("evaporative");
  if (domain_name == my_key_first) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // - pull Keys from plist
  // dependency: litter_water_content
  wc_key_ = plist_.get<std::string>("litter water content key",
          domain_name+std::string("water_content"));
  dependencies_.insert(wc_key_);

  // dependency: surface_molar_density_liquid
  rho_key_ = plist_.get<std::string>("molar density liquid key",
          domain_name+std::string("molar_density_liquid"));

  dependencies_.insert(rho_key_);

  // dependency: litter_thickness
  L_key_ = plist_.get<std::string>("litter thickness key",
          domain_name+std::string("thickness"));

  dependencies_.insert(L_key_);
}


void
EvaporativeFluxRelaxationEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> wc = S->GetPtr<CompositeVector>("litter_water_content");
  Teuchos::RCP<const CompositeVector> rho = S->GetPtr<CompositeVector>("surface_molar_density_liquid");
  Teuchos::RCP<const CompositeVector> L = S->GetPtr<CompositeVector>("litter_thickness");
  Teuchos::RCP<const CompositeVector> cv = S->GetPtr<CompositeVector>("surface_cell_volume");

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp, false);
    const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
    const Epetra_MultiVector& L_v = *L->ViewComponent(*comp, false);
    const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->EvaporativeFlux(wc_v[0][i] / cv_v[0][i], rho_v[0][i], L_v[0][i]);
    }
  }
}


void
EvaporativeFluxRelaxationEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> wc = S->GetPtr<CompositeVector>("litter_water_content");
  Teuchos::RCP<const CompositeVector> rho = S->GetPtr<CompositeVector>("surface_molar_density_liquid");
  Teuchos::RCP<const CompositeVector> L = S->GetPtr<CompositeVector>("litter_thickness");
  Teuchos::RCP<const CompositeVector> cv = S->GetPtr<CompositeVector>("surface_cell_volume");

  if (wrt_key == "litter_water_content") {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& L_v = *L->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEvaporativeFluxDLitterWaterContent(wc_v[0][i] / cv_v[0][i], rho_v[0][i], L_v[0][i]) / cv_v[0][i];
      }
    }

  } else if (wrt_key == "surface_molar_density_liquid") {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& L_v = *L->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEvaporativeFluxDSurfaceMolarDensityLiquid(wc_v[0][i] / cv_v[0][i], rho_v[0][i], L_v[0][i]);
      }
    }

  } else if (wrt_key == "litter_thickness") {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& wc_v = *wc->ViewComponent(*comp, false);
      const Epetra_MultiVector& rho_v = *rho->ViewComponent(*comp, false);
      const Epetra_MultiVector& L_v = *L->ViewComponent(*comp, false);
      const Epetra_MultiVector& cv_v = *cv->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DEvaporativeFluxDLitterThickness(wc_v[0][i] / cv_v[0][i], rho_v[0][i], L_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
