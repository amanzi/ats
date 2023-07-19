/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Satish Karra (satkarra@lanl.gov)
*/

/*
  Interface for a thermal conductivity model with three phases.

*/

#include "dbc.hh"
#include "thermal_conductivity_threephase_factory.hh"
#include "thermal_conductivity_threephase_evaluator.hh"

namespace Amanzi {
namespace Energy {

ThermalConductivityThreePhaseEvaluator::ThermalConductivityThreePhaseEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  poro_key_ = Keys::readKey(plist_, domain, "porosity", "porosity");
  dependencies_.insert(KeyTag{ poro_key_, tag });

  temp_key_ = Keys::readKey(plist_, domain, "temperature", "temperature");
  dependencies_.insert(KeyTag{ temp_key_, tag });

  sat_key_ = Keys::readKey(plist_, domain, "saturation liquid", "saturation_liquid");
  dependencies_.insert(KeyTag{ sat_key_, tag });

  sat2_key_ = Keys::readKey(plist_, domain, "second saturation key", "saturation_ice");
  dependencies_.insert(KeyTag{ sat2_key_, tag });

  AMANZI_ASSERT(plist_.isSublist("thermal conductivity parameters"));
  Teuchos::ParameterList tc_sublist = plist_.sublist("thermal conductivity parameters");

  ThermalConductivityThreePhaseFactory fac;

  for (Teuchos::ParameterList::ConstIterator lcv = tc_sublist.begin(); lcv != tc_sublist.end();
       ++lcv) {
    std::string name = lcv->first;
    if (tc_sublist.isSublist(name)) {
      Teuchos::ParameterList& tcp_sublist = tc_sublist.sublist(name);
      std::string region_name = tcp_sublist.get<std::string>("region");
      Teuchos::RCP<ThermalConductivityThreePhase> tc =
        fac.createThermalConductivityModel(tcp_sublist);
      tcs_.push_back(std::make_pair(region_name, tc));
    } else {
      Errors::Message message("ThermalConductivityThreePhaseEvaluator: region-based lists.  "
                              "(Perhaps you have an old-style input file?)");
      Exceptions::amanzi_throw(message);
    }
  }
}


Teuchos::RCP<Evaluator>
ThermalConductivityThreePhaseEvaluator::Clone() const
{
  return Teuchos::rcp(new ThermalConductivityThreePhaseEvaluator(*this));
}


void
ThermalConductivityThreePhaseEvaluator::Evaluate_(const State& S,
                                                  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> poro = S.GetPtr<CompositeVector>(poro_key_, tag);
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);
  Teuchos::RCP<const CompositeVector> sat = S.GetPtr<CompositeVector>(sat_key_, tag);
  Teuchos::RCP<const CompositeVector> sat2 = S.GetPtr<CompositeVector>(sat2_key_, tag);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result[0]->Mesh();

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    AMANZI_ASSERT(*comp == "cell");
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp, false);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp, false);
    const Epetra_MultiVector& sat2_v = *sat2->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    for (std::vector<RegionModelPair>::const_iterator lcv = tcs_.begin(); lcv != tcs_.end();
         ++lcv) {
      std::string region_name = lcv->first;
      if (mesh->isValidSetName(region_name, AmanziMesh::Entity_kind::CELL)) {
        // get the indices of the domain.
        auto id_list = mesh->getSetEntities(
          region_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

        // loop over indices
        for (const auto& id: id_list) {
          result_v[0][id] = lcv->second->ThermalConductivity(
            poro_v[0][id], sat_v[0][id], sat2_v[0][id], temp_v[0][id]);
        }
      } else {
        std::stringstream m;
        m << "Thermal conductivity evaluator: unknown region on cells: \"" << region_name << "\"";
        Errors::Message message(m.str());
        Exceptions::amanzi_throw(message);
      }
    }
  }
  result[0]->Scale(1.e-6); // convert to MJ
}


void
ThermalConductivityThreePhaseEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // pull out the dependencies
  Teuchos::RCP<const CompositeVector> poro = S.GetPtr<CompositeVector>(poro_key_, tag);
  Teuchos::RCP<const CompositeVector> temp = S.GetPtr<CompositeVector>(temp_key_, tag);
  Teuchos::RCP<const CompositeVector> sat = S.GetPtr<CompositeVector>(sat_key_, tag);
  Teuchos::RCP<const CompositeVector> sat2 = S.GetPtr<CompositeVector>(sat2_key_, tag);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = result[0]->Mesh();

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    AMANZI_ASSERT(*comp == "cell");
    const Epetra_MultiVector& poro_v = *poro->ViewComponent(*comp, false);
    const Epetra_MultiVector& temp_v = *temp->ViewComponent(*comp, false);
    const Epetra_MultiVector& sat_v = *sat->ViewComponent(*comp, false);
    const Epetra_MultiVector& sat2_v = *sat2->ViewComponent(*comp, false);
    Epetra_MultiVector& result_v = *result[0]->ViewComponent(*comp, false);

    if (wrt_key == poro_key_) {
      for (std::vector<RegionModelPair>::const_iterator lcv = tcs_.begin(); lcv != tcs_.end();
           ++lcv) {
        std::string region_name = lcv->first;
        if (mesh->isValidSetName(region_name, AmanziMesh::Entity_kind::CELL)) {
          // get the indices of the domain.
          auto id_list = mesh->getSetEntities(
            region_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

          // loop over indices
          for (const auto& id: id_list) {
            result_v[0][id] = lcv->second->DThermalConductivity_DPorosity(
              poro_v[0][id], sat_v[0][id], sat2_v[0][id], temp_v[0][id]);
          }
        } else {
          std::stringstream m;
          m << "Thermal conductivity evaluator: unknown region on cells: \"" << region_name << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }
      }

    } else if (wrt_key == sat_key_) {
      for (std::vector<RegionModelPair>::const_iterator lcv = tcs_.begin(); lcv != tcs_.end();
           ++lcv) {
        std::string region_name = lcv->first;
        if (mesh->isValidSetName(region_name, AmanziMesh::Entity_kind::CELL)) {
          // get the indices of the domain.
          auto id_list = mesh->getSetEntities(
            region_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

          // loop over indices
          for (const auto& id : id_list) {
            result_v[0][id] = lcv->second->DThermalConductivity_DSaturationLiquid(
              poro_v[0][id], sat_v[0][id], sat2_v[0][id], temp_v[0][id]);
          }
        } else {
          std::stringstream m;
          m << "Thermal conductivity evaluator: unknown region on cells: \"" << region_name << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }
      }

    } else if (wrt_key == sat2_key_) {
      for (std::vector<RegionModelPair>::const_iterator lcv = tcs_.begin(); lcv != tcs_.end();
           ++lcv) {
        std::string region_name = lcv->first;
        if (mesh->isValidSetName(region_name, AmanziMesh::Entity_kind::CELL)) {
          // get the indices of the domain.
          auto id_list = mesh->getSetEntities(
            region_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

          // loop over indices
          for (const auto& id: id_list) {
            result_v[0][id] = lcv->second->DThermalConductivity_DSaturationIce(
              poro_v[0][id], sat_v[0][id], sat2_v[0][id], temp_v[0][id]);
          }
        } else {
          std::stringstream m;
          m << "Thermal conductivity evaluator: unknown region on cells: \"" << region_name << "\"";
          Errors::Message message(m.str());
          Exceptions::amanzi_throw(message);
        }
      }

    } else if (wrt_key == temp_key_) {
      for (std::vector<RegionModelPair>::const_iterator lcv = tcs_.begin(); lcv != tcs_.end();
           ++lcv) {
        std::string region_name = lcv->first;
        if (mesh->isValidSetName(region_name, AmanziMesh::Entity_kind::CELL)) {
          // get the indices of the domain.
          auto id_list = mesh->getSetEntities(
            region_name, AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

          // loop over indices
          for (const auto& id: id_list) {
            result_v[0][id] = lcv->second->DThermalConductivity_DTemperature(
              poro_v[0][id], sat_v[0][id], sat2_v[0][id], temp_v[0][id]);
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

  result[0]->Scale(1.e-6); // convert to MJ

  Epetra_MultiVector& result_v = *result[0]->ViewComponent("cell", false);
}

} // namespace Energy
} // namespace Amanzi
