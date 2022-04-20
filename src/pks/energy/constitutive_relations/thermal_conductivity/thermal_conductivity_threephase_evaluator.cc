/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Interface for a thermal conductivity model with three phases.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
  Satish Karra (satkarra@lanl.gov)
*/

#include "dbc.hh"
#include "thermal_conductivity_threephase_factory.hh"
#include "thermal_conductivity_threephase_evaluator.hh"

namespace Amanzi {
namespace Energy {

ThermalConductivityThreePhaseEvaluator::ThermalConductivityThreePhaseEvaluator(
    Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist) {
  Key domain = Keys::getDomain(my_key_);

  poro_key_ = Keys::readKey(plist_, domain, "porosity", "porosity");
  dependencies_.insert(poro_key_);

  temp_key_ = Keys::readKey(plist_, domain, "temperature", "temperature");
  dependencies_.insert(temp_key_);

  sat_key_ = Keys::readKey(plist_, domain, "saturation liquid", "saturation_liquid");
  dependencies_.insert(sat_key_);

  sat2_key_ = Keys::readKey(plist_, domain, "second saturation key", "saturation_ice");
  dependencies_.insert(sat2_key_);

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList tc_sublist = plist_.sublist("thermal conductivity parameters");

  ThermalConductivityThreePhaseFactory fac;

  for (Teuchos::ParameterList::ConstIterator lcv=tc_sublist.begin();
       lcv!=tc_sublist.end(); ++lcv) {
    std::string name = lcv->first;
    if (tc_sublist.isSublist(name)) {
      Teuchos::ParameterList& tcp_sublist = tc_sublist.sublist(name);
      std::string region_name = tcp_sublist.get<std::string>("region");
      Teuchos::RCP<ThermalConductivityThreePhase> tc = fac.createThermalConductivityModel(tcp_sublist);
      tcs_.push_back(std::make_pair(region_name,tc));
    } else {
      Errors::Message message("ThermalConductivityThreePhaseEvaluator: region-based lists.  (Perhaps you have an old-style input file?)");
      Exceptions::amanzi_throw(message);
    }
  }
}


ThermalConductivityThreePhaseEvaluator::ThermalConductivityThreePhaseEvaluator(
    const ThermalConductivityThreePhaseEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    poro_key_(other.poro_key_),
    temp_key_(other.temp_key_),
    sat_key_(other.sat_key_),
    sat2_key_(other.sat2_key_),
    tcs_(other.tcs_) {}

Teuchos::RCP<FieldEvaluator>
ThermalConductivityThreePhaseEvaluator::Clone() const {
  return Teuchos::rcp(new ThermalConductivityThreePhaseEvaluator(*this));
}


void ThermalConductivityThreePhaseEvaluator::EvaluateField_(
    const Teuchos::Ptr<State>& S,
    const Teuchos::Ptr<CompositeVector>& result) {
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(poro_key_);
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);
  Teuchos::RCP<const CompositeVector> sat2 = S->GetFieldData(sat2_key_);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();

  for (CompositeVector::name_iterator comp = result->begin();
       comp!=result->end(); ++comp) {
    AMANZI_ASSERT(*comp == "cell");
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp,false);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp,false);
    const Epetra_MultiVector& sat2_v = *sat2->ViewComponent(*comp,false);
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
          result_v[0][*id] = lcv->second->ThermalConductivity(poro_v[0][*id],
                  sat_v[0][*id], sat2_v[0][*id], temp_v[0][*id]);
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


void ThermalConductivityThreePhaseEvaluator::EvaluateFieldPartialDerivative_(
    const Teuchos::Ptr<State>& S, Key wrt_key,
    const Teuchos::Ptr<CompositeVector>& result) {
  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> poro = S->GetFieldData(poro_key_);
  Teuchos::RCP<const CompositeVector> temp = S->GetFieldData(temp_key_);
  Teuchos::RCP<const CompositeVector> sat = S->GetFieldData(sat_key_);
  Teuchos::RCP<const CompositeVector> sat2 = S->GetFieldData(sat2_key_);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();

  for (CompositeVector::name_iterator comp = result->begin();
       comp!=result->end(); ++comp) {
    AMANZI_ASSERT(*comp == "cell");
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp,false);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp,false);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp,false);
    const Epetra_MultiVector& sat2_v = *sat2->ViewComponent(*comp,false);
    Epetra_MultiVector& result_v = *result->ViewComponent(*comp,false);

    if (wrt_key == poro_key_) {
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
            result_v[0][*id] = lcv->second->DThermalConductivity_DPorosity(poro_v[0][*id],
                                              sat_v[0][*id], sat2_v[0][*id], temp_v[0][*id]);
          }
        } else {
          std::stringstream m;
          m << "Thermal conductivity evaluator: unknown region on cells: \"" << region_name << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }
      }
      
    } else if (wrt_key == sat_key_) {
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
            result_v[0][*id] = lcv->second->DThermalConductivity_DSaturationLiquid(
                                    poro_v[0][*id], sat_v[0][*id], sat2_v[0][*id], temp_v[0][*id]);
          }
        } else {
          std::stringstream m;
          m << "Thermal conductivity evaluator: unknown region on cells: \"" << region_name << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }
      }

    } else if (wrt_key == sat2_key_) {
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
            result_v[0][*id] = lcv->second->DThermalConductivity_DSaturationIce(
                                    poro_v[0][*id], sat_v[0][*id], sat2_v[0][*id], temp_v[0][*id]);
          }
        } else {
          std::stringstream m;
          m << "Thermal conductivity evaluator: unknown region on cells: \"" << region_name << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }
      }

    } else if (wrt_key == temp_key_) {
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
            result_v[0][*id] = lcv->second->DThermalConductivity_DTemperature(
                                    poro_v[0][*id], sat_v[0][*id], sat2_v[0][*id], temp_v[0][*id]);
          }
        } else {
          std::stringstream m;
          m << "Thermal conductivity evaluator: unknown region on cells: \"" << region_name << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }
      }
      
    } else {
      AMANZI_ASSERT(false);
    }
  }
    
  result->Scale(1.e-6); // convert to MJ

  Epetra_MultiVector& result_v = *result->ViewComponent("cell",false);
}


} //namespace
} //namespace

