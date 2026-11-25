/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Downregulates bare soil evaporation via the Lee-Pielke cosine ramp on
  liquid saturation relative to a WRM-derived field capacity.
*/

#include <cmath>

#include "MeshAlgorithms.hh"
#include "evaporation_downregulation_lee_pielke_evaluator.hh"

namespace Amanzi {
namespace Flow {

EvaporationDownregulationLeePielkeEvaluator::EvaporationDownregulationLeePielkeEvaluator(
  Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  Key domain_surf = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;
  Key domain_ss = Keys::readDomainHint(plist, domain_surf, "surface", "subsurface");

  pot_evap_key_ =
    Keys::readKey(plist_, domain_surf, "potential evaporation", "potential_evaporation");
  dependencies_.insert(KeyTag{ pot_evap_key_, tag });

  sat_gas_key_ = Keys::readKey(plist_, domain_ss, "gas saturation", "saturation_gas");
  dependencies_.insert(KeyTag{ sat_gas_key_, tag });

  std::string params_name = plist_.get<std::string>("model parameters", "WRM parameters");
  Teuchos::ParameterList& sublist = plist_.sublist(params_name);
  wrms_ = createWRMPartition(sublist);
}

Teuchos::RCP<Evaluator>
EvaporationDownregulationLeePielkeEvaluator::Clone() const
{
  return Teuchos::rcp(new EvaporationDownregulationLeePielkeEvaluator(*this));
}


void
EvaporationDownregulationLeePielkeEvaluator::Evaluate_(const State& S,
                                                        const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result[0]->Mesh()->getParentMesh(), -1);
    wrms_->first->Verify();
  }

  Teuchos::RCP<const CompositeVector> pot_evap = S.GetPtr<CompositeVector>(pot_evap_key_, tag);
  Teuchos::RCP<const CompositeVector> sat_gas = S.GetPtr<CompositeVector>(sat_gas_key_, tag);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = sat_gas->Mesh();
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = result[0]->Mesh();

  for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end(); ++comp) {
    AMANZI_ASSERT(*comp == "cell"); // partition on cell only
    const Epetra_MultiVector& pot_evap_v = *(pot_evap->ViewComponent(*comp, false));
    const Epetra_MultiVector& sat_gas_v = *(sat_gas->ViewComponent(*comp, false));
    Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

    int count = result[0]->size(*comp);
    for (int sc = 0; sc != count; ++sc) {
      int f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);

      double sat_fc = wrms_->second[(*wrms_->first)[c]]->field_capacity();
      double one_minus_sg = 1.0 - sat_gas_v[0][c];

      double beta;
      if (one_minus_sg >= sat_fc) {
        beta = 1.0;
      } else {
        double arg = M_PI * one_minus_sg / sat_fc;
        beta = 0.25 * std::pow(1.0 - std::cos(arg), 2.0);
      }

      double e_pot = pot_evap_v[0][sc];
      result_v[0][sc] = (e_pot > 0.0) ? e_pot * beta : e_pot;
    }
  }
}


void
EvaporationDownregulationLeePielkeEvaluator::EvaluatePartialDerivative_(
  const State& S,
  const Key& wrt_key,
  const Tag& wrt_tag,
  const std::vector<CompositeVector*>& result)
{
  Tag tag = my_keys_.front().second;

  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result[0]->Mesh()->getParentMesh(), -1);
    wrms_->first->Verify();
  }

  Teuchos::RCP<const CompositeVector> pot_evap = S.GetPtr<CompositeVector>(pot_evap_key_, tag);
  Teuchos::RCP<const CompositeVector> sat_gas = S.GetPtr<CompositeVector>(sat_gas_key_, tag);
  Teuchos::RCP<const AmanziMesh::Mesh> mesh = sat_gas->Mesh();
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = result[0]->Mesh();

  if (wrt_key == pot_evap_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      AMANZI_ASSERT(*comp == "cell");
      const Epetra_MultiVector& pot_evap_v = *(pot_evap->ViewComponent(*comp, false));
      const Epetra_MultiVector& sat_gas_v = *(sat_gas->ViewComponent(*comp, false));
      Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

      int count = result[0]->size(*comp);
      for (int sc = 0; sc != count; ++sc) {
        int f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);

        double sat_fc = wrms_->second[(*wrms_->first)[c]]->field_capacity();
        double one_minus_sg = 1.0 - sat_gas_v[0][c];

        double beta;
        if (one_minus_sg >= sat_fc) {
          beta = 1.0;
        } else {
          double arg = M_PI * one_minus_sg / sat_fc;
          beta = 0.25 * std::pow(1.0 - std::cos(arg), 2.0);
        }

        double e_pot = pot_evap_v[0][sc];
        result_v[0][sc] = (e_pot > 0.0) ? beta : 1.0;
      }
    }

  } else if (wrt_key == sat_gas_key_) {
    for (CompositeVector::name_iterator comp = result[0]->begin(); comp != result[0]->end();
         ++comp) {
      AMANZI_ASSERT(*comp == "cell");
      const Epetra_MultiVector& pot_evap_v = *(pot_evap->ViewComponent(*comp, false));
      const Epetra_MultiVector& sat_gas_v = *(sat_gas->ViewComponent(*comp, false));
      Epetra_MultiVector& result_v = *(result[0]->ViewComponent(*comp, false));

      int count = result[0]->size(*comp);
      for (int sc = 0; sc != count; ++sc) {
        int f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
        int c = AmanziMesh::getFaceOnBoundaryInternalCell(*mesh, f);

        double sat_fc = wrms_->second[(*wrms_->first)[c]]->field_capacity();
        double one_minus_sg = 1.0 - sat_gas_v[0][c];
        double e_pot = pot_evap_v[0][sc];

        double dbeta_dsg;
        if (one_minus_sg >= sat_fc) {
          dbeta_dsg = 0.0;
        } else {
          // beta = 0.25 * (1 - cos(pi * (1 - sg) / sfc))^2
          // d(1-sg)/dsg = -1
          double arg = M_PI * one_minus_sg / sat_fc;
          double dbeta_darg = 0.5 * (1.0 - std::cos(arg)) * std::sin(arg);
          double darg_dsg = -M_PI / sat_fc;
          dbeta_dsg = dbeta_darg * darg_dsg;
        }

        result_v[0][sc] = (e_pot > 0.0) ? e_pot * dbeta_dsg : 0.0;
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}


void
EvaporationDownregulationLeePielkeEvaluator::EnsureCompatibility_ToDeps_(State& S)
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
