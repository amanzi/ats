/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  ViscosityEvaluator is the interface between state/data and the model, a VPM.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "pc_ice_water.hh"
#include "pc_ice_evaluator.hh"

namespace Amanzi {
namespace Flow {

PCIceEvaluator::PCIceEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {

  // my keys
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("capillary pressure of ice-water key",
            "capillary_pressure_liq_ice");
  }

  Key domain_name = Keys::getDomain(my_key_);
  // -- temperature
 
 temp_key_ = plist_.get<std::string>("temperature key", Keys::getKey(domain_name,"temperature"));

  dependencies_.insert(temp_key_);

  // Construct my PCIce model
  model_ = Teuchos::rcp(new PCIceWater(plist_.sublist("capillary pressure of ice-water")));

  if (model_->IsMolarBasis()) {
    dens_key_ = plist_.get<std::string>("molar density key", Keys::getKey(domain_name,"molar_density_liquid"));
  } else {
    dens_key_ = plist_.get<std::string>("mass density key", Keys::getKey(domain_name, "mass_density_liquid"));

  }
  dependencies_.insert(dens_key_);
};


PCIceEvaluator::PCIceEvaluator(const PCIceEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    model_(other.model_),
    temp_key_(other.temp_key_),
    dens_key_(other.dens_key_) {}


Teuchos::RCP<FieldEvaluator> PCIceEvaluator::Clone() const {
  return Teuchos::rcp(new PCIceEvaluator(*this));
}


void PCIceEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
                         const Teuchos::Ptr<CompositeVector>& result) {
  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);
  Teuchos::RCP<const CompositeVector> dens = S->GetPtr<CompositeVector>(dens_key_);
  double lambda = S->HasField("continuation_parameter") ?
    std::pow(10., -2*(*S->GetScalarData("continuation_parameter"))) : 1.;
  

  // evaluate pc
  for (CompositeVector::name_iterator comp=result->begin();
       comp!=result->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
    const Epetra_MultiVector& dens_v = *(dens->ViewComponent(*comp,false));
    Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

    int count = result->size(*comp);
    for (int id=0; id!=count; ++id) {
      result_v[0][id] = lambda * model_->CapillaryPressure(temp_v[0][id], dens_v[0][id]);
    }
  }
}


void PCIceEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {

  // Pull dependencies out of state.
  Teuchos::RCP<const CompositeVector> temp = S->GetPtr<CompositeVector>(temp_key_);
  Teuchos::RCP<const CompositeVector> dens = S->GetPtr<CompositeVector>(dens_key_);
  double lambda = S->HasField("continuation_parameter") ?
    std::pow(10., -2*(*S->GetScalarData("continuation_parameter"))) : 1.;

  if (wrt_key == temp_key_) {
    // evaluate d/dT( pc )
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
      const Epetra_MultiVector& dens_v = *(dens->ViewComponent(*comp,false));
      Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

      int count = result->size(*comp);
      for (int id=0; id!=count; ++id) {
        result_v[0][id] = lambda * model_->DCapillaryPressureDT(temp_v[0][id], dens_v[0][id]);
      }
    }
  } else if (wrt_key == dens_key_) {
    // evaluate d/drho( pc )
    for (CompositeVector::name_iterator comp=result->begin();
         comp!=result->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *(temp->ViewComponent(*comp,false));
      const Epetra_MultiVector& dens_v = *(dens->ViewComponent(*comp,false));
      Epetra_MultiVector& result_v = *(result->ViewComponent(*comp,false));

      int count = result->size(*comp);
      for (int id=0; id!=count; ++id) {
        result_v[0][id] = lambda * model_->DCapillaryPressureDRho(temp_v[0][id], dens_v[0][id]);
      }
    }
  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace
} // namespace
