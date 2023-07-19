/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  This WRM evaluator evaluates saturation of gas, liquid, and ice from
  capillary pressures for the ice-liquid and liquid-gas pairs.

*/

#include "wrm_permafrost_evaluator.hh"
#include "wrm_partition.hh"

namespace Amanzi {
namespace Flow {

/* --------------------------------------------------------------------------------
  Constructor from just a ParameterList, reads WRMs and permafrost models from list.
 -------------------------------------------------------------------------------- */
WRMPermafrostEvaluator::WRMPermafrostEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist)
{
  // get the WRMs
  AMANZI_ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList wrm_plist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(wrm_plist);

  // and the permafrost models
  AMANZI_ASSERT(plist_.isSublist("permafrost model parameters"));
  Teuchos::ParameterList perm_plist = plist_.sublist("permafrost model parameters");
  permafrost_models_ = createWRMPermafrostModelPartition(perm_plist, wrms_);

  InitializeFromPlist_();
}


/* --------------------------------------------------------------------------------
  Constructor with WRMs.
 -------------------------------------------------------------------------------- */
WRMPermafrostEvaluator::WRMPermafrostEvaluator(Teuchos::ParameterList& plist,
                                               const Teuchos::RCP<WRMPartition>& wrms)
  : EvaluatorSecondaryMonotypeCV(plist), wrms_(wrms)
{
  // and the permafrost models
  AMANZI_ASSERT(plist_.isSublist("permafrost model parameters"));
  Teuchos::ParameterList perm_plist = plist_.sublist("permafrost model parameters");
  permafrost_models_ = createWRMPermafrostModelPartition(perm_plist, wrms_);

  InitializeFromPlist_();
}


/* --------------------------------------------------------------------------------
  Constructor with Permafrost models.
 -------------------------------------------------------------------------------- */
WRMPermafrostEvaluator::WRMPermafrostEvaluator(
  Teuchos::ParameterList& plist,
  const Teuchos::RCP<WRMPermafrostModelPartition>& models)
  : EvaluatorSecondaryMonotypeCV(plist), permafrost_models_(models)
{
  InitializeFromPlist_();
}


/* --------------------------------------------------------------------------------
  Virtual opy constructor as a Evaluator.
 -------------------------------------------------------------------------------- */
Teuchos::RCP<Evaluator>
WRMPermafrostEvaluator::Clone() const
{
  return Teuchos::rcp(new WRMPermafrostEvaluator(*this));
}


/* --------------------------------------------------------------------------------
  Initialization of keys.
 -------------------------------------------------------------------------------- */
void
WRMPermafrostEvaluator::InitializeFromPlist_()
{
  // my keys are for saturation -- order matters... gas -> liq -> ice
  Key akey = my_keys_.front().first;
  Key domain_name = Keys::getDomain(akey);
  akey = Keys::getVarName(akey);
  Tag tag = my_keys_.front().second;
  my_keys_.clear();

  Key gaskey, liqkey, icekey;
  std::size_t liq_pos = akey.find("liquid");
  std::size_t ice_pos = akey.find("ice");
  std::size_t gas_pos = akey.find("gas");
  if (liq_pos != std::string::npos) {
    liqkey = Keys::readKey(plist_, domain_name, "liquid saturation", akey);
    gaskey = akey.substr(0, liq_pos) + "gas" + akey.substr(liq_pos + 6);
    gaskey = Keys::readKey(plist_, domain_name, "gas saturation", gaskey);
    icekey = akey.substr(0, liq_pos) + "ice" + akey.substr(liq_pos + 6);
    icekey = Keys::readKey(plist_, domain_name, "ice saturation", icekey);

  } else if (ice_pos != std::string::npos) {
    icekey = Keys::readKey(plist_, domain_name, "ice saturation", akey);
    gaskey = akey.substr(0, ice_pos) + "gas" + akey.substr(ice_pos + 3);
    gaskey = Keys::readKey(plist_, domain_name, "gas saturation", gaskey);
    liqkey = akey.substr(0, ice_pos) + "liquid" + akey.substr(ice_pos + 3);
    liqkey = Keys::readKey(plist_, domain_name, "liquid saturation", liqkey);

  } else if (gas_pos != std::string::npos) {
    gaskey = Keys::readKey(plist_, domain_name, "gas saturation", akey);
    icekey = akey.substr(0, gas_pos) + "ice" + akey.substr(gas_pos + 3);
    icekey = Keys::readKey(plist_, domain_name, "ice saturation", icekey);
    liqkey = akey.substr(0, gas_pos) + "liquid" + akey.substr(gas_pos + 3);
    liqkey = Keys::readKey(plist_, domain_name, "liquid saturation", liqkey);

  } else {
    liqkey = Keys::readKey(plist_, domain_name, "liquid saturation", "saturation_liquid");
    icekey = Keys::readKey(plist_, domain_name, "ice saturation", "saturation_ice");
    gaskey = Keys::readKey(plist_, domain_name, "gas saturation", "saturation_gas");
  }

  my_keys_.emplace_back(KeyTag{ gaskey, tag });
  my_keys_.emplace_back(KeyTag{ liqkey, tag });
  my_keys_.emplace_back(KeyTag{ icekey, tag });

  // liquid-gas capillary pressure
  pc_liq_key_ = Keys::readKey(
    plist_, domain_name, "gas-liquid capillary pressure", "capillary_pressure_gas_liq");
  dependencies_.insert(KeyTag{ pc_liq_key_, tag });

  // liquid-ice capillary pressure
  pc_ice_key_ = Keys::readKey(
    plist_, domain_name, "liquid-ice capillary pressure", "capillary_pressure_liq_ice");
  dependencies_.insert(KeyTag{ pc_ice_key_, tag });
}


void
WRMPermafrostEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  // Initialize the MeshPartition
  if (!permafrost_models_->first->initialized()) {
    permafrost_models_->first->Initialize(results[0]->Mesh(), -1);
    permafrost_models_->first->Verify();
  }

  // Cell values
  Epetra_MultiVector& satg_c = *results[0]->ViewComponent("cell", false);
  Epetra_MultiVector& satl_c = *results[1]->ViewComponent("cell", false);
  Epetra_MultiVector& sati_c = *results[2]->ViewComponent("cell", false);

  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& pc_liq_c =
    *S.GetPtr<CompositeVector>(pc_liq_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& pc_ice_c =
    *S.GetPtr<CompositeVector>(pc_ice_key_, tag)->ViewComponent("cell", false);

  double sats[3];
  int ncells = satg_c.MyLength();
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    int i = (*permafrost_models_->first)[c];
    permafrost_models_->second[i]->saturations(pc_liq_c[0][c], pc_ice_c[0][c], sats);
    satg_c[0][c] = sats[0];
    satl_c[0][c] = sats[1];
    sati_c[0][c] = sats[2];
  }

  // Potentially do face values as well, though only for saturation_liquid?
  if (results[0]->HasComponent("boundary_face")) {
    Epetra_MultiVector& satg_bf = *results[0]->ViewComponent("boundary_face", false);
    Epetra_MultiVector& satl_bf = *results[1]->ViewComponent("boundary_face", false);
    Epetra_MultiVector& sati_bf = *results[2]->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pc_liq_bf =
      *S.GetPtr<CompositeVector>(pc_liq_key_, tag)->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pc_ice_bf =
      *S.GetPtr<CompositeVector>(pc_ice_key_, tag)->ViewComponent("boundary_face", false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = results[0]->Mesh();
    const Epetra_Map& vandelay_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);
    const Epetra_Map& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);

    // calculate boundary face values
    int nbfaces = satg_bf.MyLength();
    for (int bf = 0; bf != nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);

      int i = (*permafrost_models_->first)[cells[0]];
      permafrost_models_->second[i]->saturations(pc_liq_bf[0][bf], pc_ice_bf[0][bf], sats);
      satg_bf[0][bf] = sats[0];
      satl_bf[0][bf] = sats[1];
      sati_bf[0][bf] = sats[2];
    }
  }
}


void
WRMPermafrostEvaluator::EvaluatePartialDerivative_(const State& S,
                                                   const Key& wrt_key,
                                                   const Tag& wrt_tag,
                                                   const std::vector<CompositeVector*>& results)
{
  // Initialize the MeshPartition
  if (!permafrost_models_->first->initialized()) {
    permafrost_models_->first->Initialize(results[0]->Mesh(), -1);
    permafrost_models_->first->Verify();
  }

  // Cell values
  Epetra_MultiVector& satg_c = *results[0]->ViewComponent("cell", false);
  Epetra_MultiVector& satl_c = *results[1]->ViewComponent("cell", false);
  Epetra_MultiVector& sati_c = *results[2]->ViewComponent("cell", false);

  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& pc_liq_c =
    *S.GetPtr<CompositeVector>(pc_liq_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& pc_ice_c =
    *S.GetPtr<CompositeVector>(pc_ice_key_, tag)->ViewComponent("cell", false);

  double dsats[3];
  if (wrt_key == pc_liq_key_) {
    int ncells = satg_c.MyLength();
    for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
      int i = (*permafrost_models_->first)[c];
      permafrost_models_->second[i]->dsaturations_dpc_liq(pc_liq_c[0][c], pc_ice_c[0][c], dsats);

      satg_c[0][c] = dsats[0];
      satl_c[0][c] = dsats[1];
      sati_c[0][c] = dsats[2];
    }

  } else if (wrt_key == pc_ice_key_) {
    int ncells = satg_c.MyLength();
    for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
      int i = (*permafrost_models_->first)[c];
      permafrost_models_->second[i]->dsaturations_dpc_ice(pc_liq_c[0][c], pc_ice_c[0][c], dsats);

      satg_c[0][c] = dsats[0];
      satl_c[0][c] = dsats[1];
      sati_c[0][c] = dsats[2];
    }
  } else {
    AMANZI_ASSERT(0);
  }

  // Potentially do face values as well, though only for saturation_liquid?
  if (results[0]->HasComponent("boundary_face")) {
    Epetra_MultiVector& satg_bf = *results[0]->ViewComponent("boundary_face", false);
    Epetra_MultiVector& satl_bf = *results[1]->ViewComponent("boundary_face", false);
    Epetra_MultiVector& sati_bf = *results[2]->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pc_liq_bf =
      *S.GetPtr<CompositeVector>(pc_liq_key_, tag)->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pc_ice_bf =
      *S.GetPtr<CompositeVector>(pc_ice_key_, tag)->ViewComponent("boundary_face", false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = results[0]->Mesh();
    const Epetra_Map& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);
    const Epetra_Map& vandelay_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);

    if (wrt_key == pc_liq_key_) {
      // calculate boundary face values
      int nbfaces = satl_bf.MyLength();
      for (int bf = 0; bf != nbfaces; ++bf) {
        // given a boundary face, we need the internal cell to choose the right WRM
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        AMANZI_ASSERT(cells.size() == 1);

        int i = (*permafrost_models_->first)[cells[0]];
        permafrost_models_->second[i]->dsaturations_dpc_liq(
          pc_liq_bf[0][bf], pc_ice_bf[0][bf], dsats);
        satg_bf[0][bf] = dsats[0];
        satl_bf[0][bf] = dsats[1];
        sati_bf[0][bf] = dsats[2];
      }

    } else if (wrt_key == pc_ice_key_) {
      // calculate boundary face values
      int nbfaces = satl_bf.MyLength();
      for (int bf = 0; bf != nbfaces; ++bf) {
        // given a boundary face, we need the internal cell to choose the right WRM
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
        AMANZI_ASSERT(cells.size() == 1);

        int i = (*permafrost_models_->first)[cells[0]];
        permafrost_models_->second[i]->dsaturations_dpc_ice(
          pc_liq_bf[0][bf], pc_ice_bf[0][bf], dsats);
        satg_bf[0][bf] = dsats[0];
        satl_bf[0][bf] = dsats[1];
        sati_bf[0][bf] = dsats[2];
      }
    } else {
      AMANZI_ASSERT(0);
    }
  }
}


} // namespace Flow
} // namespace Amanzi
