/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  ViscosityEvaluator is the interface between state/data and the model, a VPM.

*/

#include "viscosity_relation_factory.hh"
#include "viscosity_evaluator.hh"

namespace Amanzi {
namespace Relations {

ViscosityEvaluator::ViscosityEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  // Set up my dependencies.
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- temperature
  temp_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });

  // Construct my Viscosity model
  AMANZI_ASSERT(plist_.isSublist("viscosity model parameters"));
  ViscosityRelationFactory visc_fac;
  visc_ = visc_fac.createViscosity(plist_.sublist("viscosity model parameters"));
};


Teuchos::RCP<Evaluator>
ViscosityEvaluator::Clone() const
{
  return Teuchos::rcp(new ViscosityEvaluator(*this));
}


void
ViscosityEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  // Pull dependencies out of state.
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);

  // evaluate p_s / p_atm
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int id = 0; id != count; ++id) {
      AMANZI_ASSERT(temp_v[0][id] > 200.);
      result_v[0][id] = visc_->Viscosity(temp_v[0][id]);
    }
  }
}


void
ViscosityEvaluator::EvaluatePartialDerivative_(const State& S,
                                               const Key& wrt_key,
                                               const Tag& wrt_tag,
                                               const std::vector<CompositeVector*>& result)
{
  AMANZI_ASSERT(wrt_key == temp_key_);

  // Pull dependencies out of state.
  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);

  // evaluate d/dT( p_s / p_atm )
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int id = 0; id != count; ++id) { result_v[0][id] = visc_->DViscosityDT(temp_v[0][id]); }
  }
}

} // namespace Relations
} // namespace Amanzi
