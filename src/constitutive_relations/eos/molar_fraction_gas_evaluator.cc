/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Determining the molar fraction of a gas component within a gas mixture.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "vapor_pressure_relation_factory.hh"
#include "molar_fraction_gas_evaluator.hh"

namespace Amanzi {
namespace Relations {

MolarFractionGasEvaluator::MolarFractionGasEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  // set up the actual model
  AMANZI_ASSERT(plist_.isSublist("vapor pressure model parameters"));
  VaporPressureRelationFactory vpm_fac;
  sat_vapor_model_ = vpm_fac.createVaporPressure(
      plist_.sublist("vapor pressure model parameters"));

  // process the list for my provided field.
  if (my_key_ == "")
    my_key_ = plist_.get<std::string>("molar fraction key");

  // set up dependencies
  std::size_t end = my_key_.find_first_of("_");
  std::string domain_name = my_key_.substr(0,end);
  if (domain_name == std::string("mol") ||
      domain_name == std::string("molar")) {
    domain_name = std::string("");
  } else {
    domain_name = domain_name+std::string("_");
  }
  domain_name =Keys::getDomain(my_key_);
  
  temp_key_= plist_.get<std::string>("temperature key",
                                     Keys::getKey(domain_name,"temperature"));
  dependencies_.insert(temp_key_);
}


MolarFractionGasEvaluator::MolarFractionGasEvaluator(const MolarFractionGasEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    sat_vapor_model_(other.sat_vapor_model_),
    temp_key_(other.temp_key_) {}


Teuchos::RCP<FieldEvaluator>
MolarFractionGasEvaluator::Clone() const {
  return Teuchos::rcp(new MolarFractionGasEvaluator(*this));
}


void MolarFractionGasEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));

  // evaluate p_s / p_atm
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      AMANZI_ASSERT(temp_v[0][id] > 200.);
      result_v[0][id] = sat_vapor_model_->SaturatedVaporPressure(temp_v[0][id]) / p_atm;
    }
  }
}


void MolarFractionGasEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  AMANZI_ASSERT(wrt_key == temp_key_);

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);
  const double& p_atm = *(S->GetScalarData("atmospheric_pressure"));

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = sat_vapor_model_->DSaturatedVaporPressureDT(temp_v[0][id]) / p_atm;
    }
  }
}


} // namespace
} // namespace

