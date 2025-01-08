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

#include "pc_ice_water.hh"
#include "pc_ice_evaluator.hh"

namespace Amanzi {
namespace Flow {

PCIceEvaluator::PCIceEvaluator(Teuchos::ParameterList& plist) : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  temp_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });

  // Construct my PCIce model
  model_ = Teuchos::rcp(new PCIceWater(plist_.sublist("capillary pressure of ice-water")));
  if (model_->IsMolarBasis()) {
    dens_key_ = Keys::readKey(plist_, domain_name, "molar density", "molar_density_liquid");
    dependencies_.insert(KeyTag{ dens_key_, tag });
  } else {
    dens_key_ = Keys::readKey(plist_, domain_name, "mass density", "mass_density_liquid");
    dependencies_.insert(KeyTag{ dens_key_, tag });
  }
};


Teuchos::RCP<Evaluator>
PCIceEvaluator::Clone() const
{
  return Teuchos::rcp(new PCIceEvaluator(*this));
}


void
PCIceEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);
  Teuchos::RCP<const CompositeVector> dens = S.GetPtr<CompositeVector>(dens_key_, tag);
  double lambda = S.HasData<double>("continuation_parameter", Tags::DEFAULT) ?
                    std::pow(10., -2 * (S.Get<double>("continuation_parameter", Tags::DEFAULT))) :
                    1.;

  // evaluate pc
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp, false));
    const Epetra_MultiVector& dens_v = *(dens->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int id = 0; id != count; ++id) {
      result_v[0][id] = lambda * model_->CapillaryPressure(temp_v[0][id], dens_v[0][id]);
    }
  }
}


void
PCIceEvaluator::EvaluatePartialDerivative_(const State& S,
                                           const Key& wrt_key,
                                           const Tag& wrt_tag,
                                           const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);
  Teuchos::RCP<const CompositeVector> dens = S.GetPtr<CompositeVector>(dens_key_, tag);
  double lambda = S.HasData<double>("continuation_parameter", Tags::DEFAULT) ?
                    std::pow(10., -2 * (S.Get<double>("continuation_parameter", Tags::DEFAULT))) :
                    1.;

  if (wrt_key == temp_key_) {
    // evaluate d/dT( pc )
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp, false));
      const Epetra_MultiVector& dens_v = *(dens->ViewComponent(*comp, false));
      Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

      int count = result[0]->size(*comp);
      for (int id = 0; id != count; ++id) {
        result_v[0][id] = lambda * model_->DCapillaryPressureDT(temp_v[0][id], dens_v[0][id]);
      }
    }
  } else if (wrt_key == dens_key_) {
    // evaluate d/drho( pc )
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp, false));
      const Epetra_MultiVector& dens_v = *(dens->ViewComponent(*comp, false));
      Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

      int count = result[0]->size(*comp);
      for (int id = 0; id != count; ++id) {
        result_v[0][id] = lambda * model_->DCapillaryPressureDRho(temp_v[0][id], dens_v[0][id]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace Flow
} // namespace Amanzi
