/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/* -----------------------------------------------------------------------------
ATS

Authors: Svetlana Tokareva (tokareva@lanl.gov)

FieldEvaluator for water density.
----------------------------------------------------------------------------- */


#include "lake_energy_evaluator.hh"

namespace Amanzi {
namespace LakeThermo {

LakeEnergyEvaluator::LakeEnergyEvaluator(Teuchos::ParameterList& plist) :
    EvaluatorSecondaryMonotypeCV(plist) {

  // Set up my dependencies.
  std::string domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- temperature
  temperature_key_ = Keys::readKey(plist_, domain_name, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temperature_key_, tag });

  // -- density
  density_key_ = Keys::readKey(plist_, domain_name, "density", "density");
  dependencies_.insert(KeyTag{ density_key_, tag });

  // -- heat capacity
  heat_capacity_key_ = Keys::readKey(plist_, domain_name, "heat capacity", "heat_capacity");
  dependencies_.insert(KeyTag{ heat_capacity_key_, tag });

};

Teuchos::RCP<Evaluator>
LakeEnergyEvaluator::Clone() const {
  return Teuchos::rcp(new LakeEnergyEvaluator(*this));
};


void LakeEnergyEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{

  Tag tag = my_keys_.front().second;

  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

  // evaluate density
  const Epetra_MultiVector& rho_v =
      *S.GetPtr<CompositeVector>(density_key_,tag)->ViewComponent("cell",false);

  // evaluate heat capacity
  const Epetra_MultiVector& cp_v =
      *S.GetPtr<CompositeVector>(heat_capacity_key_,tag)->ViewComponent("cell",false);

  for (CompositeVector::name_iterator comp=result[0]->begin();
       comp!=result[0]->end(); ++comp) {
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

    int ncomp = result[0]->size(*comp, false);
    for (int i=0; i!=ncomp; ++i) {
      double T = temp_v[0][i];
      double rho = rho_v[0][i];
      double cp = cp_v[0][i];
      double L_f = -333500.; // plus or minus??? // !!!!!add this is cp*rho*T terms in advection fluxzes and sources!!!!
      if (T < 273.15) {
        result_v[0][i] = rho*(cp*T + 1.*L_f);
      } else {
        result_v[0][i] = rho*(cp*T - 0.*L_f);
      }

    }
  }
};


void LakeEnergyEvaluator::EvaluatePartialDerivative_(const State& S,
                                              const Key& wrt_key,
                                              const Tag& wrt_tag,
                                              const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  result[0]->PutScalar(0.);

  if (wrt_key == temperature_key_) {
    Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temperature_key_,tag);

    // evaluate density
    const Epetra_MultiVector& rho_v =
        *S.GetPtr<CompositeVector>(density_key_,tag)->ViewComponent("cell",false);

    // evaluate heat capacity
    const Epetra_MultiVector& cp_v =
        *S.GetPtr<CompositeVector>(heat_capacity_key_,tag)->ViewComponent("cell",false);

    for (CompositeVector::name_iterator comp=result[0]->begin();
         comp!=result[0]->end(); ++comp) {
      const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
      Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp,false);

      int ncomp = result[0]->size(*comp, false);
      for (int i=0; i!=ncomp; ++i) {
        double T = temp_v[0][i];
        double rho = rho_v[0][i];
        double cp = cp_v[0][i];
        result_v[0][i] = rho*cp;
      }
    }
  }
};


} //namespace
} //namespace
