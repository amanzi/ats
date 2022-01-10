/*
  The latent heat from evaporative flux evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    keyDeclarationList =   Key qe_key_;
    myKeyFirst = latent
    evalNameString = latent heat from evaporative flux
    keyCopyConstructorList =     qe_key_(other.qe_key_),
    evalNameCaps = LATENT_HEAT
    namespaceCaps = SURFACEBALANCE
    paramDeclarationList =   double Le_;
    modelDerivDeclarationList =   double DLatentHeatDEvaporativeFlux(double qe) const;
    evalClassName = LatentHeat
    keyCompositeVectorList =   Teuchos::RCP<const CompositeVector> qe = S->GetPtr<CompositeVector>("evaporative_flux");
    namespace = SurfaceBalance
    modelInitializeParamsList =   Le_ = plist.get<double>("latent heat of vaporization [MJ/mol]", 0.0449994810744);
    myMethodDeclarationArgs = double qe
    myKey = latent_heat
    evalName = latent_heat
    modelMethodDeclaration =   double LatentHeat(double qe) const;
    myKeyMethod = LatentHeat
    myMethodArgs = qe_v[0][i]  
  
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "latent_heat_evaluator.hh"
#include "latent_heat_model.hh"

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

// Constructor from ParameterList
LatentHeatEvaluator::LatentHeatEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("latent_heat parameters");
  model_ = Teuchos::rcp(new LatentHeatModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
LatentHeatEvaluator::LatentHeatEvaluator(const LatentHeatEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    qe_key_(other.qe_key_),    
    model_(other.model_) {}


// Virtual copy constructor
Teuchos::RCP<FieldEvaluator>
LatentHeatEvaluator::Clone() const
{
  return Teuchos::rcp(new LatentHeatEvaluator(*this));
}


// Initialize by setting up dependencies
void
LatentHeatEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);

  std::string my_key_first("latent");
  if (domain_name == my_key_first) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }

  // - pull Keys from plist
  // dependency: evaporative_flux
  qe_key_ = plist_.get<std::string>("evaporative flux key",
          domain_name+std::string("evaporative_flux"));
  dependencies_.insert(qe_key_);
}


void
LatentHeatEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> qe = S->GetPtr<CompositeVector>("evaporative_flux");

  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& qe_v = *qe->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    int ncomp = result->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      result_v[0][i] = model_->LatentHeat(qe_v[0][i]);
    }
  }
}


void
LatentHeatEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> qe = S->GetPtr<CompositeVector>("evaporative_flux");

  if (wrt_key == "evaporative_flux") {
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& qe_v = *qe->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

      int ncomp = result->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        result_v[0][i] = model_->DLatentHeatDEvaporativeFlux(qe_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} //namespace
} //namespace
} //namespace
