/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ViscosityEvaluator is the interface between state/data and the model, a VPM.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "viscosity_relation_factory.hh"
#include "viscosity_evaluator.hh"

namespace Amanzi {
namespace Relations {

ViscosityEvaluator::ViscosityEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondary(plist) {

  // my keys
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("viscosity key", "viscosity_liquid");
  }

  // Set up my dependencies.
  Key domain_name = Keys::getDomain(my_key_);

  // -- temperature
  temp_key_ = plist_.get<std::string>("temperature key",
          Keys::getKey(domain_name, "temperature"));
  dependencies_.emplace_back(std::make_pair(temp_key_, my_tag_));

  // Construct my Viscosity model
  ASSERT(plist_.isSublist("viscosity model parameters"));
  ViscosityRelationFactory visc_fac;
  visc_ = visc_fac.createViscosity(plist_.sublist("viscosity model parameters"));
};


ViscosityEvaluator::ViscosityEvaluator(const ViscosityEvaluator& other) :
    EvaluatorSecondary(other),
    visc_(other.visc_),
    temp_key_(other.temp_key_) {}


Teuchos::RCP<Evaluator> ViscosityEvaluator::Clone() const {
  return Teuchos::rcp(new ViscosityEvaluator(*this));
}


void ViscosityEvaluator::Evaluate_(const State& S,
                                   CompositeVector& result) {
  // Pull dependencies out of state.
  const CompositeVector& temp = S.Get<CompositeVector>(temp_key_);

  // evaluate p_s / p_atm
  for (CompositeVector::name_iterator comp=result.begin(); comp!=result.end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp.ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result.ViewComponent(*comp,false));

    int count = result.size(*comp);
    for (int id=0; id!=count; ++id) {
      ASSERT(temp_v[0][id] > 200.);
      result_v[0][id] = visc_->Viscosity(temp_v[0][id]);
    }
  }
}


void ViscosityEvaluator::EvaluatePartialDerivative_(
                                                    const State& S, const Key& wrt_key, const Key& wrt_tag,
                                                    CompositeVector& result) {
  ASSERT(wrt_key == temp_key_);
  ASSERT(wrt_tag == my_tag_);
  
  // Pull dependencies out of state.
  const CompositeVector& temp = S.Get<CompositeVector>(temp_key_);

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp=result.begin();
       comp!=result.end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp.ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result.ViewComponent(*comp,false));

    int count = result.size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = visc_->DViscosityDT(temp_v[0][id]);
    }
  }
}

} // namespace
} // namespace
