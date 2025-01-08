/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  The ideal gas equation of state evaluator is an algebraic evaluator of a given model.

  Generated via evaluator_generator with:
    modelInitializeParamsList =   cv_ = plist.get<double>("heat capacity");
    evalName = eos_ideal_gas
    modelMethodDeclaration =   double Density(double temp, double pres) const;
    namespaceCaps = GENERAL
    namespace = General
    paramDeclarationList =   double cv_;
    evalNameCaps = EOS_IDEAL_GAS
    myMethodArgs = temp_v[0][i], pres_v[0][i]
    myKeyMethod = Density
    myKeyFirst = density
    evalNameString = ideal gas equation of state
    myMethodDeclarationArgs = double temp, double pres
    evalClassName = EosIdealGas
    myKey = density

*/

#include "eos_ideal_gas_evaluator.hh"
#include "eos_ideal_gas_model.hh"

namespace Amanzi {
namespace General {
namespace Relations {

// Constructor from ParameterList
EosIdealGasEvaluator::EosIdealGasEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Teuchos::ParameterList& sublist = plist_.sublist("eos_ideal_gas parameters");
  model_ = Teuchos::rcp(new EosIdealGasModel(sublist));
  InitializeFromPlist_();
}


// Copy constructor
EosIdealGasEvaluator::EosIdealGasEvaluator(const EosIdealGasEvaluator& other)
  : EvaluatorSecondaryMonotypeCV(other),
    temp_key_(other.temp_key_),
    pres_key_(other.pres_key_),
    model_(other.model_)
{}


// Virtual copy constructor
Teuchos::RCP<Evaluator>
EosIdealGasEvaluator::Clone() const
{
  return Teuchos::rcp(new EosIdealGasEvaluator(*this));
}


// Initialize by setting up dependencies
void
EosIdealGasEvaluator::InitializeFromPlist_()
{
  // Set up my dependencies
  // - defaults to prefixed via domain
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0, end);

  std::string my_key_first("density");
  if (domain_name == my_key_first) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name + std::string("_");
  }

  // - pull Keys from plist
  // dependency: temperature
  temp_key_ = plist_.get<std::string>("temperature key", domain_name + std::string("temperature"));
  dependencies_.insert(temp_key_);

  // dependency: pressure
  pres_key_ = plist_.get<std::string>("pressure key", domain_name + std::string("pressure"));
  dependencies_.insert(pres_key_);
}


void
EosIdealGasEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                                     const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

    int ncomp = result->size(*comp, false);
    for (int i = 0; i != ncomp; ++i) {
      result_v[0][i] = model_->Density(temp_v[0][i], pres_v[0][i]);
    }
  }
}


void
EosIdealGasEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
                                                      Key wrt_key,
                                                      const Teuchos::Ptr<CompositeVector>& result)
{
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> pres = S->GetFieldData(pres_key_);

  if (wrt_key == temp_key_) {
    for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

      int ncomp = result->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DDensityDTemperature(temp_v[0][i], pres_v[0][i]);
      }
    }

  } else if (wrt_key == pres_key_) {
    for (CompositeVector::name_iterator comp = result->begin(); comp != result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
      const Epetra_MultiVector& pres_v = *pres->ViewComponent(*comp, false);
      Epetra_MultiVector& result_v = *result->ViewComponent(*comp, false);

      int ncomp = result->size(*comp, false);
      for (int i = 0; i != ncomp; ++i) {
        result_v[0][i] = model_->DDensityDPressure(temp_v[0][i], pres_v[0][i]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


} // namespace Relations
} // namespace General
} // namespace Amanzi
