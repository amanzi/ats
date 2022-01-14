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
    EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  temp_key_= Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{temp_key_, tag});

  // set up the actual model
  VaporPressureRelationFactory vpm_fac;
  sat_vapor_model_ = vpm_fac.createVaporPressure(
      plist_.sublist("vapor pressure model parameters"));

}

Teuchos::RCP<Evaluator>
MolarFractionGasEvaluator::Clone() const {
  return Teuchos::rcp(new MolarFractionGasEvaluator(*this));
}


void MolarFractionGasEvaluator::Evaluate_(const State& S,
        const std::vector<CompositeVector*>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_);
  const double& p_atm = S.Get<double>("atmospheric_pressure");

  // evaluate p_s / p_atm
  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp,false));

    int count = result[0]->size(*comp);
    for (int id=0; id!=count; ++id) {
      AMANZI_ASSERT(temp_v[0][id] > 200.);
      result_v[0][id] = sat_vapor_model_->SaturatedVaporPressure(temp_v[0][id]) / p_atm;
    }
  }
}


void MolarFractionGasEvaluator::EvaluatePartialDerivative_(
    const State& S, const Key& wrt_key, const Tag& wrt_tag,
    const std::vector<CompositeVector*>& result) {
  AMANZI_ASSERT(wrt_key == temp_key_);

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_);
  const double& p_atm = S.Get<double>("atmospheric_pressure");

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp,false));

    int count = result[0]->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = sat_vapor_model_->DSaturatedVaporPressureDT(temp_v[0][id]) / p_atm;
    }
  }
}


} // namespace
} // namespace

