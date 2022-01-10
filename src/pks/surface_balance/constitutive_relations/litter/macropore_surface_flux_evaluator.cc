/*
  The macropore-surface flux evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    evalName = macropore_surface_flux
    modelMethodDeclaration =   double MacroporeSurfaceFlux(double pM, double ps, double krs, double krM, double K) const;
    namespaceCaps = SURFACEBALANCE
    namespace = SurfaceBalance
    evalNameCaps = MACROPORE_SURFACE_FLUX
    myMethodArgs = pM_v[0][i], ps_v[0][i], krs_v[0][i], krM_v[0][i], K_v[0][i]
    myKeyMethod = MacroporeSurfaceFlux
    myKeyFirst = macropore
    evalNameString = macropore-surface flux
    myMethodDeclarationArgs = double pM, double ps, double krs, double krM, double K
    evalClassName = MacroporeSurfaceFlux
    myKey = macropore_surface_flux  
  
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "macropore_surface_flux_evaluator.hh"
#include "macropore_surface_flux_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
MacroporeSurfaceFluxEvaluator::MacroporeSurfaceFluxEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("macropore_surface_flux parameters");
  model_ = Teuchos::rcp(new MacroporeSurfaceFluxModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
MacroporeSurfaceFluxEvaluator::MacroporeSurfaceFluxEvaluator(const MacroporeSurfaceFluxEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    pM_key_(other.pM_key_),
    ps_key_(other.ps_key_),
    krs_key_(other.krs_key_),
    krM_key_(other.krM_key_),
    K_key_(other.K_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
MacroporeSurfaceFluxEvaluator::Clone() const
{
  return Teuchos::rcp(new MacroporeSurfaceFluxEvaluator(*this));
}


// Initialize by setting up dependencies
void
MacroporeSurfaceFluxEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);

  std::string my_key_first("macropore");
  if (domain_name == my_key_first) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // - pull Keys from plist
  // dependency: macropore_pressure
  pM_key_ = plist_.get<std::string>("macropore pressure key",
          domain_name+std::string("macropore_pressure"));
  dependencies_.insert(pM_key_);

  // dependency: surface_pressure
  ps_key_ = plist_.get<std::string>("surface pressure key",
          domain_name+std::string("surface_pressure"));
  dependencies_.insert(ps_key_);

  // dependency: surface_relative_permeability
  krs_key_ = plist_.get<std::string>("surface relative permeability key",
          domain_name+std::string("surface_relative_permeability"));
  dependencies_.insert(krs_key_);

  // dependency: macropore_relative_permeability
  krM_key_ = plist_.get<std::string>("macropore relative permeability key",
          domain_name+std::string("macropore_relative_permeability"));
  dependencies_.insert(krM_key_);

  // dependency: macropore_absolute_permeability
  K_key_ = plist_.get<std::string>("macropore absolute permeability key",
          domain_name+std::string("macropore_absolute_permeability"));
  dependencies_.insert(K_key_);
}


void
MacroporeSurfaceFluxEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> pM = S->GetPtr<CompositeVector>(pM_key_);
Teuchos::RCP<const CompositeVector> ps = S->GetPtr<CompositeVector>(ps_key_);
Teuchos::RCP<const CompositeVector> krs = S->GetPtr<CompositeVector>(krs_key_);
Teuchos::RCP<const CompositeVector> krM = S->GetPtr<CompositeVector>(krM_key_);
Teuchos::RCP<const CompositeVector> K = S->GetPtr<CompositeVector>(K_key_);

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
    const Epetra_MultiVector& ps_v = *ps->ViewComponent(*comp, false);
    const Epetra_MultiVector& krs_v = *krs->ViewComponent(*comp, false);
    const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
    const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->MacroporeSurfaceFlux(pM_v[0][i], ps_v[0][i], krs_v[0][i], krM_v[0][i], K_v[0][i]);        
    }
  }
}


void
MacroporeSurfaceFluxEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
Teuchos::RCP<const CompositeVector> pM = S->GetPtr<CompositeVector>(pM_key_);
Teuchos::RCP<const CompositeVector> ps = S->GetPtr<CompositeVector>(ps_key_);
Teuchos::RCP<const CompositeVector> krs = S->GetPtr<CompositeVector>(krs_key_);
Teuchos::RCP<const CompositeVector> krM = S->GetPtr<CompositeVector>(krM_key_);
Teuchos::RCP<const CompositeVector> K = S->GetPtr<CompositeVector>(K_key_);

  if (wrt_key == pM_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& ps_v = *ps->ViewComponent(*comp, false);
      const Epetra_MultiVector& krs_v = *krs->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMacroporeSurfaceFluxDMacroporePressure(pM_v[0][i], ps_v[0][i], krs_v[0][i], krM_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == ps_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& ps_v = *ps->ViewComponent(*comp, false);
      const Epetra_MultiVector& krs_v = *krs->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMacroporeSurfaceFluxDSurfacePressure(pM_v[0][i], ps_v[0][i], krs_v[0][i], krM_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == krs_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& ps_v = *ps->ViewComponent(*comp, false);
      const Epetra_MultiVector& krs_v = *krs->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMacroporeSurfaceFluxDSurfaceRelativePermeability(pM_v[0][i], ps_v[0][i], krs_v[0][i], krM_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == krM_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& ps_v = *ps->ViewComponent(*comp, false);
      const Epetra_MultiVector& krs_v = *krs->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMacroporeSurfaceFluxDMacroporeRelativePermeability(pM_v[0][i], ps_v[0][i], krs_v[0][i], krM_v[0][i], K_v[0][i]);
      }
    }

  } else if (wrt_key == K_key_) {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& pM_v = *pM->ViewComponent(*comp, false);
      const Epetra_MultiVector& ps_v = *ps->ViewComponent(*comp, false);
      const Epetra_MultiVector& krs_v = *krs->ViewComponent(*comp, false);
      const Epetra_MultiVector& krM_v = *krM->ViewComponent(*comp, false);
      const Epetra_MultiVector& K_v = *K->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DMacroporeSurfaceFluxDMacroporeAbsolutePermeability(pM_v[0][i], ps_v[0][i], krs_v[0][i], krM_v[0][i], K_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
