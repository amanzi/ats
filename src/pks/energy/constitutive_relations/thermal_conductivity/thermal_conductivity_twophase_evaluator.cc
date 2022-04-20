/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with two phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "dbc.hh"
#include "thermal_conductivity_twophase_factory.hh"
#include "thermal_conductivity_twophase_evaluator.hh"

namespace Amanzi {
namespace Energy {

ThermalConductivityTwoPhaseEvaluator::ThermalConductivityTwoPhaseEvaluator(
      Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  Key domain = Keys::getDomain(my_key_);

  poro_key_ = Keys::readKey(plist_, domain, "porosity", "porosity");
  dependencies_.insert(poro_key_);

  sat_key_ = Keys::readKey(plist_, domain, "saturation liquid", "saturation_liquid");
  dependencies_.insert(sat_key_);


  Teuchos::ParameterList tc_sublist = plist_.sublist("thermal conductivity parameters");
  ThermalConductivityTwoPhaseFactory fac;
  for (Teuchos::ParameterList::ConstIterator lcv=tc_sublist.begin();
       lcv!=tc_sublist.end(); ++lcv) {
    std::string name = lcv->first;
    if (tc_sublist.isSublist(name)) {
      Teuchos::ParameterList& tcp_sublist = tc_sublist.sublist(name);
      std::string region_name = tcp_sublist.get<std::string>("region");
      Teuchos::RCP<ThermalConductivityTwoPhase> tc = fac.createThermalConductivityModel(tcp_sublist);
      tcs_.push_back(std::make_pair(region_name,tc));
    } else {
      Errors::Message message("ThermalConductivityTwoPhaseEvaluator: region-based lists.  (Perhaps you have an old-style input file?)");
      Exceptions::amanzi_throw(message);
    }
  }
}


ThermalConductivityTwoPhaseEvaluator::ThermalConductivityTwoPhaseEvaluator(
      const ThermalConductivityTwoPhaseEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    poro_key_(other.poro_key_),
    sat_key_(other.sat_key_),
    tcs_(other.tcs_) {}

Teuchos::RCP<FieldEvaluator>
ThermalConductivityTwoPhaseEvaluator::Clone() const {
  return Teuchos::rcp(new ThermalConductivityTwoPhaseEvaluator(*this));
}


void ThermalConductivityTwoPhaseEvaluator::EvaluateField_(
      const Teuchos::Ptr<State>& S,
      const Teuchos::Ptr<CompositeVector>& result) {
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(poro_key_);
  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();

  for (CompositeVector::name_iterator comp = result->begin();
       comp!=result->end(); ++comp) {
    AMANZI_ASSERT(*comp == "cell");
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp,false);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    for (std::vector<RegionModelPair>::const_iterator lcv = tcs_.begin();
         lcv != tcs_.end(); ++lcv) {
      std::string region_name = lcv->first;
      if (mesh->valid_set_name(region_name, AmanziMesh::CELL)) {
        // get the indices of the domain.
        AmanziMesh::Entity_ID_List id_list;
        mesh->get_set_entities(region_name, AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED, &id_list);

        // loop over indices
        for (AmanziMesh::Entity_ID_List::const_iterator id=id_list.begin();
             id!=id_list.end(); ++id) {
          result_v[0][*id] = lcv->second->ThermalConductivity(poro_v[0][*id], sat_v[0][*id]);
        }
      } else {
        std::stringstream m;
        m << "Thermal conductivity evaluator: unknown region on cells: \"" << region_name << "\"";
        Errors::Message message(m.str());
        Exceptions::amanzi_throw(message);
      }
    }
  }
  result->Scale(1.e-6); // convert to MJ
}


void ThermalConductivityTwoPhaseEvaluator::EvaluateFieldPartialDerivative_(
      const Teuchos::Ptr<State>& S, Key wrt_key,
      const Teuchos::Ptr<CompositeVector>& result) {
  // not yet implemented in underlying models!
  result->PutScalar(0.);
  // result->Scale(1.e-6); // convert to MJ
}

} //namespace
} //namespace
