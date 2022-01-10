/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  The WRM Evaluator simply calls the WRM with the correct arguments.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


#include "wrm_evaluator.hh"
#include "wrm_factory.hh"

namespace Amanzi {
namespace Flow {

WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariablesFieldEvaluator(plist),
    calc_other_sat_(true) {

  AMANZI_ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList wrm_plist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(wrm_plist);

  InitializeFromPlist_();
}

WRMEvaluator::WRMEvaluator(Teuchos::ParameterList& plist,
                           const Teuchos::RCP<WRMPartition>& wrms) :
    SecondaryVariablesFieldEvaluator(plist),
    wrms_(wrms) {
  InitializeFromPlist_();
}

WRMEvaluator::WRMEvaluator(const WRMEvaluator& other) :
    SecondaryVariablesFieldEvaluator(other),
    calc_other_sat_(other.calc_other_sat_),
    cap_pres_key_(other.cap_pres_key_),
    wrms_(other.wrms_) {}


Teuchos::RCP<FieldEvaluator> WRMEvaluator::Clone() const {
  return Teuchos::rcp(new WRMEvaluator(*this));
}

void WRMEvaluator::InitializeFromPlist_() {
  // my keys are for saturation, note that order matters, liquid -> gas
  Key akey = Keys::cleanPListName(plist_.name());
  Key domain_name = Keys::getDomain(akey);

  std::size_t liq_pos = akey.find("liquid");
  if (liq_pos != std::string::npos) {
    my_keys_.push_back(plist_.get<std::string>("saturation key", akey));
    Key otherkey = akey.substr(0,liq_pos)+"gas"+akey.substr(liq_pos+6);
    my_keys_.push_back(plist_.get<std::string>("other saturation key", otherkey));
  } else {
    std::size_t gas_pos = akey.find("gas");
    if (gas_pos != std::string::npos) {
      Key otherkey = akey.substr(0,gas_pos)+"liquid"+akey.substr(gas_pos+3);
      my_keys_.push_back(plist_.get<std::string>("saturation key", otherkey));
      my_keys_.push_back(plist_.get<std::string>("other saturation key", akey));
    } else {
      my_keys_.push_back(Keys::readKey(plist_, domain_name,
              "saturation", "saturation_liquid"));
      my_keys_.push_back(Keys::readKey(plist_, domain_name,
              "other saturation", "saturation_gas"));
    }
  }

  // my dependencies are capillary pressure.
  cap_pres_key_ = Keys::readKey(plist_, domain_name, "capillary pressure key",
          "capillary_pressure_gas_liq");
  dependencies_.insert(cap_pres_key_);
}


void WRMEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const std::vector<Teuchos::Ptr<CompositeVector> >& results) {
  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(results[0]->Mesh(), -1);
    wrms_->first->Verify();
  }

  Epetra_MultiVector& sat_c = *results[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *S->GetPtr<CompositeVector>(cap_pres_key_)
      ->ViewComponent("cell",false);

  // calculate cell values
  AmanziMesh::Entity_ID ncells = sat_c.MyLength();
  for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
    sat_c[0][c] = wrms_->second[(*wrms_->first)[c]]->saturation(pres_c[0][c]);
  }

  // Potentially do face values as well.
  if (results[0]->HasComponent("boundary_face")) {
    Epetra_MultiVector& sat_bf = *results[0]->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& pres_bf = *S->GetPtr<CompositeVector>(cap_pres_key_)
        ->ViewComponent("boundary_face",false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = results[0]->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
    const Epetra_Map& face_map = mesh->face_map(false);
    AmanziMesh::Entity_ID_List cells;

    // calculate boundary face values
    int nbfaces = sat_bf.MyLength();
    for (int bf=0; bf!=nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      int index = (*wrms_->first)[cells[0]];
      sat_bf[0][bf] = wrms_->second[index]->saturation(pres_bf[0][bf]);
    }
  }

  // If needed, also do gas saturation
  if (calc_other_sat_) {
    for (CompositeVector::name_iterator comp=results[1]->begin();
         comp!=results[1]->end(); ++comp) {

      if (results[0]->HasComponent(*comp)) {
        // sat_g = 1 - sat_l
        results[1]->ViewComponent(*comp,false)->PutScalar(1.);
        results[1]->ViewComponent(*comp,false)->Update(-1,
                *results[0]->ViewComponent(*comp,false), 1.);
      } else {
        // sat_l not available on this component, loop and call the model

        // Currently this is not ever the case.  If this error shows up, it
        // can easily be implemented. -- etc
        AMANZI_ASSERT(0);
      }
    }
  }
}


void WRMEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const std::vector<Teuchos::Ptr<CompositeVector> > & results) {
  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(results[0]->Mesh(), -1);
    wrms_->first->Verify();
  }

  AMANZI_ASSERT(wrt_key == cap_pres_key_);

  Epetra_MultiVector& sat_c = *results[0]->ViewComponent("cell",false);
  const Epetra_MultiVector& pres_c = *S->GetPtr<CompositeVector>(cap_pres_key_)
      ->ViewComponent("cell",false);

  // calculate cell values
  AmanziMesh::Entity_ID ncells = sat_c.MyLength();
  for (AmanziMesh::Entity_ID c=0; c!=ncells; ++c) {
    sat_c[0][c] = wrms_->second[(*wrms_->first)[c]]->d_saturation(pres_c[0][c]);
  }

  // Potentially do face values as well.
  if (results[0]->HasComponent("boundary_face")) {
    Epetra_MultiVector& sat_bf = *results[0]->ViewComponent("boundary_face",false);
    const Epetra_MultiVector& pres_bf = *S->GetPtr<CompositeVector>(cap_pres_key_)
        ->ViewComponent("boundary_face",false);

    // Need to get boundary face's inner cell to specify the WRM.
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = results[0]->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
    const Epetra_Map& face_map = mesh->face_map(false);
    AmanziMesh::Entity_ID_List cells;

    // calculate boundary face values
    int nbfaces = sat_bf.MyLength();
    for (int bf=0; bf!=nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      int index = (*wrms_->first)[cells[0]];
      sat_bf[0][bf] = wrms_->second[index]->d_saturation(pres_bf[0][bf]);
    }
  }

  // If needed, also do gas saturation
  if (calc_other_sat_) {
    for (CompositeVector::name_iterator comp=results[1]->begin();
         comp!=results[1]->end(); ++comp) {

      if (results[0]->HasComponent(*comp)) {
        // d_sat_g =  - d_sat_l
        results[1]->ViewComponent(*comp,false)->Update(-1,
                *results[0]->ViewComponent(*comp,false), 0.);
      } else {
        // sat_l not available on this component, loop and call the model

        // Currently this is not ever the case.  If this error shows up, it
        // can easily be implemented. -- etc
        AMANZI_ASSERT(0);
      }
    }
  }

}


} //namespace
} //namespace
