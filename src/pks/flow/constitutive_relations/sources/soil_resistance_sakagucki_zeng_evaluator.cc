/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Bo Gao (gaob@ornl.gov)
*/

/*
  Evaluates the soil resistance at top cells through the Sakagucki-Zeng
  method and assign them to surface cells.
*/


#include "MeshAlgorithms.hh"
#include "soil_resistance_sakagucki_zeng_evaluator.hh"
#include "soil_resistance_sakagucki_zeng_model.hh"

namespace Amanzi {
namespace Flow {

SoilResistanceSakaguckiZengEvaluator::SoilResistanceSakaguckiZengEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  std::string domain_surf = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  Key domain_ss = Keys::readDomainHint(plist, domain_surf, "surface", "subsurface");

  sat_gas_key_ = Keys::readKey(plist_, domain_ss, "gas saturation", "saturation_gas");
  dependencies_.insert(KeyTag{ sat_gas_key_, tag });

  poro_key_ = Keys::readKey(plist_, domain_ss, "porosity", "porosity");
  dependencies_.insert(KeyTag{ poro_key_, tag });

  std::string params_name = plist_.get<std::string>("model parameters", "WRM parameters");
  Teuchos::ParameterList& sublist = plist_.sublist(params_name);
  models_ = createSoilResistanceModelPartition(sublist);
}

Teuchos::RCP<Evaluator>
SoilResistanceSakaguckiZengEvaluator::Clone() const
{
  return Teuchos::rcp(new SoilResistanceSakaguckiZengEvaluator(*this));
}


// Required methods from EvaluatorSecondaryMonotypeCV
void
SoilResistanceSakaguckiZengEvaluator::Evaluate_(const State& S,
                                                const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;
  // Initialize the MeshPartition
  if (!models_->first->initialized()) {
    models_->first->Initialize(result[0]->Mesh()->getParentMesh(), -1);
    models_->first->Verify();
  }

  Teuchos::RCP<const CompositeVector> sat_gas = S.GetPtr<CompositeVector>(sat_gas_key_, tag);
  Teuchos::RCP<const CompositeVector> poro = S.GetPtr<CompositeVector>(poro_key_, tag);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = sat_gas->Mesh();
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = result[0]->Mesh();

  // evaluate the model
  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    AMANZI_ASSERT(*comp == "cell"); // partition on cell only
    const Epetra_MultiVector& sat_gas_v = *(sat_gas->ViewComponent(*comp, false));
    const Epetra_MultiVector& poro_v = *(poro->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int sc = 0; sc != count; ++sc) {
      int f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
      result_v[0][sc] =
        models_->second[(*models_->first)[c]]->RsoilbySakagickiZeng(sat_gas_v[0][c], poro_v[0][c]);
    }
  }
}


void
SoilResistanceSakaguckiZengEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  // Initialize the MeshPartition
  if (!models_->first->initialized()) {
    models_->first->Initialize(result[0]->Mesh(), -1);
    models_->first->Verify();
  }

  Tag tag = my_keys_.front().second;
  Teuchos::RCP<const CompositeVector> sat_gas = S.GetPtr<CompositeVector>(sat_gas_key_, tag);
  Teuchos::RCP<const CompositeVector> poro = S.GetPtr<CompositeVector>(poro_key_, tag);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = sat_gas->Mesh();
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = result[0]->Mesh();

  if (wrt_key == sat_gas_key_) {
    // evaluate the model
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      AMANZI_ASSERT(*comp == "cell"); // partition on cell only
      const Epetra_MultiVector& sat_gas_v = *(sat_gas->ViewComponent(*comp, false));
      const Epetra_MultiVector& poro_v = *(poro->ViewComponent(*comp, false));
      Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

      int count = result[0]->size(*comp);
      for (int sc = 0; sc != count; ++sc) {
        int f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        result_v[0][sc] = models_->second[(*models_->first)[c]]->DRsoilbySakagickiZengDSatGas(
          sat_gas_v[0][c], poro_v[0][c]);
      }
    }

  } else if (wrt_key == poro_key_) {
    // evaluate the model
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      AMANZI_ASSERT(*comp == "cell"); // partition on cell only
      const Epetra_MultiVector& sat_gas_v = *(sat_gas->ViewComponent(*comp, false));
      const Epetra_MultiVector& poro_v = *(poro->ViewComponent(*comp, false));
      Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

      int count = result[0]->size(*comp);
      for (int sc = 0; sc != count; ++sc) {
        int f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);
        result_v[0][sc] = models_->second[(*models_->first)[c]]->DRsoilbySakagickiZengDPorosity(
          sat_gas_v[0][c], poro_v[0][c]);
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


void
SoilResistanceSakaguckiZengEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  const auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
                                                                     my_keys_.front().second);
  if (fac.Mesh() != Teuchos::null) {
    CompositeVectorSpace dep_fac;
    dep_fac.SetMesh(fac.Mesh()->getParentMesh())
      ->SetGhosted(true)
      ->AddComponent("cell", AmanziMesh::CELL, 1);

    for (const auto& dep : dependencies_) {
      if (Keys::getDomain(dep.first) == Keys::getDomain(my_keys_.front().first)) {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(fac);
      } else {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(dep_fac);
      }
    }
  }
}


} // namespace Flow
} // namespace Amanzi
