/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

*/


#include "wrm_evaluator.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {

WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), calc_other_sat_(true)
{
  AMANZI_ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList wrm_plist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(wrm_plist);

  InitializeFromPlist_();
}

WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist, const Teuchos::RCP<WRMPartition>& wrms)
  : EvaluatorSecondaryMonotypeCV(plist), wrms_(wrms)
{
  InitializeFromPlist_();
}


Teuchos::RCP<Evaluator>
WRMEvaluator::Clone() const
{
  return Teuchos::rcp(new WRMEvaluator(*this));
}

void
WRMEvaluator::InitializeFromPlist_()
{
  // my keys are for saturation, note that order matters, liquid -> gas
  Key akey = my_keys_.front().first;
  Key domain_name = Keys::getDomain(akey);
  Tag tag = my_keys_.front().second;
  my_keys_.clear();

  std::size_t liq_pos = akey.find("liquid");
  std::size_t gas_pos = akey.find("gas");
  if (liq_pos != std::string::npos) {
    my_keys_.emplace_back(KeyTag{ akey, tag });

    Key otherkey = akey.substr(0, liq_pos) + "gas" + akey.substr(liq_pos + 6);
    otherkey = Keys::readKey(plist_, domain_name, "other saturation", otherkey);
    my_keys_.emplace_back(KeyTag{ otherkey, tag });

  } else if (gas_pos != std::string::npos) {
    Key otherkey = akey.substr(0, gas_pos) + "liquid" + akey.substr(gas_pos + 3);
    otherkey = Keys::readKey(plist_, domain_name, "saturation", otherkey);
    my_keys_.emplace_back(KeyTag{ otherkey, tag });
    my_keys_.emplace_back(KeyTag{ akey, tag });

  } else {
    Key liquid_key = Keys::readKey(plist_, domain_name, "saturation");
    Key gas_key = Keys::readKey(plist_, domain_name, "other saturation");
    my_keys_.emplace_back(KeyTag{ liquid_key, tag });
    my_keys_.emplace_back(KeyTag{ gas_key, tag });
  }

  // my dependencies are capillary pressure.
  cap_pres_key_ =
    Keys::readKey(plist_, domain_name, "capillary pressure key", "capillary_pressure_gas_liq");
  dependencies_.insert(KeyTag{ cap_pres_key_, tag });
}


void
WRMEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(results[0]->Mesh(), -1);
    wrms_->first->Verify();
  }

  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& sat_c = *results[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c =
    *S.GetPtr<CompositeVector>(cap_pres_key_, tag)->ViewComponent("cell", false);

  // calculate cell values
  AmanziMesh::Entity_ID ncells = sat_c.MyLength();
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    sat_c[0][c] = wrms_->second[(*wrms_->first)[c]]->saturation(pres_c[0][c]);
  }

  // Potentially do face values as well.
  if (results[0]->HasComponent("boundary_face")) {
    Epetra_MultiVector& sat_bf = *results[0]->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pres_bf =
      *S.GetPtr<CompositeVector>(cap_pres_key_, tag)->ViewComponent("boundary_face", false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = results[0]->Mesh();
    const Epetra_Map& vandelay_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);
    const Epetra_Map& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);

    // calculate boundary face values
    int nbfaces = sat_bf.MyLength();
    for (int bf = 0; bf != nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);

      int index = (*wrms_->first)[cells[0]];
      sat_bf[0][bf] = wrms_->second[index]->saturation(pres_bf[0][bf]);
    }
  }

  // If needed, also do gas saturation
  if (calc_other_sat_) {
    for (CompositeVector::name_iterator comp = results[1]->begin(); comp != results[1]->end();
         ++comp) {
      if (results[0]->HasComponent(*comp)) {
        // sat_g = 1 - sat_l
        results[1]->ViewComponent(*comp, false)->PutScalar(1.);
        results[1]
          ->ViewComponent(*comp, false)
          ->Update(-1, *results[0]->ViewComponent(*comp, false), 1.);
      } else {
        // sat_l not available on this component, loop and call the model

        // Currently this is not ever the case.  If this error shows up, it
        // can easily be implemented. -- etc
        AMANZI_ASSERT(0);
      }
    }
  }
}


void
WRMEvaluator::EvaluatePartialDerivative_(const State& S,
                                         const Key& wrt_key,
                                         const Tag& wrt_tag,
                                         const std::vector<CompositeVector*>& results)
{
  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(results[0]->Mesh(), -1);
    wrms_->first->Verify();
  }

  Tag tag = my_keys_.front().second;
  Epetra_MultiVector& sat_c = *results[0]->ViewComponent("cell", false);
  const Epetra_MultiVector& pres_c =
    *S.GetPtr<CompositeVector>(cap_pres_key_, tag)->ViewComponent("cell", false);

  // calculate cell values
  AmanziMesh::Entity_ID ncells = sat_c.MyLength();
  for (AmanziMesh::Entity_ID c = 0; c != ncells; ++c) {
    sat_c[0][c] = wrms_->second[(*wrms_->first)[c]]->d_saturation(pres_c[0][c]);
  }

  // Potentially do face values as well.
  if (results[0]->HasComponent("boundary_face")) {
    Epetra_MultiVector& sat_bf = *results[0]->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& pres_bf =
      *S.GetPtr<CompositeVector>(cap_pres_key_, tag)->ViewComponent("boundary_face", false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = results[0]->Mesh();
    const Epetra_Map& vandelay_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);
    const Epetra_Map& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);

    // calculate boundary face values
    int nbfaces = sat_bf.MyLength();
    for (int bf = 0; bf != nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);

      int index = (*wrms_->first)[cells[0]];
      sat_bf[0][bf] = wrms_->second[index]->d_saturation(pres_bf[0][bf]);
    }
  }

  // If needed, also do gas saturation
  if (calc_other_sat_) {
    for (CompositeVector::name_iterator comp = results[1]->begin(); comp != results[1]->end();
         ++comp) {
      if (results[0]->HasComponent(*comp)) {
        // d_sat_g =  - d_sat_l
        results[1]
          ->ViewComponent(*comp, false)
          ->Update(-1, *results[0]->ViewComponent(*comp, false), 0.);
      } else {
        // sat_l not available on this component, loop and call the model

        // Currently this is not ever the case.  If this error shows up, it
        // can easily be implemented. -- etc
        AMANZI_ASSERT(0);
      }
    }
  }
}


} // namespace Flow
} // namespace Amanzi
