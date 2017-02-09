/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  Rel perm( pc ( sat ) ).

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "rel_perm_evaluator.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    min_val_(0.) {

  ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(sublist);
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMPartition>& wrms) :
    SecondaryVariableFieldEvaluator(plist),
    wrms_(wrms),
    min_val_(0.) {
  InitializeFromPlist_();
}

RelPermEvaluator::RelPermEvaluator(const RelPermEvaluator& other) :
    SecondaryVariableFieldEvaluator(other),
    wrms_(other.wrms_),
    sat_key_(other.sat_key_),
    dens_key_(other.dens_key_),
    visc_key_(other.visc_key_),
    surf_rel_perm_key_(other.surf_rel_perm_key_),
    is_dens_visc_(other.is_dens_visc_),
    is_surf_(other.is_surf_),
    surf_mesh_key_(other.surf_mesh_key_),
    perm_scale_(other.perm_scale_),
    min_val_(other.min_val_) {}


Teuchos::RCP<FieldEvaluator>
RelPermEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermEvaluator(*this));
}


void RelPermEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("rel perm key", "relative_permeability");
  }

  // dependencies
  Key domain_name = getDomain(my_key_);
  
  // -- saturation liquid
  sat_key_ = plist_.get<std::string>("saturation key",
          getKey(domain_name, "saturation_liquid"));
  dependencies_.insert(sat_key_);

  is_dens_visc_ = plist_.get<bool>("use density on viscosity in rel perm", true);
  if (is_dens_visc_) {
    dens_key_ = plist_.get<std::string>("density key",
            getKey(domain_name, "molar_density_liquid"));
    dependencies_.insert(dens_key_);

    visc_key_ = plist_.get<std::string>("viscosity key",
            getKey(domain_name, "viscosity_liquid"));
    dependencies_.insert(visc_key_);
  }

  // surface alterations
  is_surf_ = plist_.get<bool>("use surface rel perm", false);
  if (is_surf_) {
    surf_rel_perm_key_ = plist_.get<std::string>("surface rel perm key",
            "surface-relative_permeability");
    dependencies_.insert(surf_rel_perm_key_);

    surf_mesh_key_ = plist_.get<std::string>("surface mesh key", "surface");
  }
  
  // cutoff above 0?
  min_val_ = plist_.get<double>("minimum rel perm cutoff", 0.);
  perm_scale_ = plist_.get<double>("permeability rescaling");
}


// Special purpose EnsureCompatibility required because of surface rel perm.
void RelPermEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  if (!is_surf_) {
    SecondaryVariableFieldEvaluator::EnsureCompatibility(S);
  } else {

    // Ensure my field exists.  Requirements should be already set.
    ASSERT(my_key_ != std::string(""));
    Teuchos::RCP<CompositeVectorSpace> my_fac = S->RequireField(my_key_, my_key_);

    // check plist for vis or checkpointing control
    bool io_my_key = plist_.get<bool>(std::string("visualize ")+my_key_, true);
    S->GetField(my_key_, my_key_)->set_io_vis(io_my_key);
    bool checkpoint_my_key = plist_.get<bool>(std::string("checkpoint ")+my_key_, false);
    S->GetField(my_key_, my_key_)->set_io_checkpoint(checkpoint_my_key);

    // If my requirements have not yet been set, we'll have to hope they
    // get set by someone later.  For now just defer.
    if (my_fac->Mesh() != Teuchos::null) {
      // Create an unowned factory to check my dependencies.
      Teuchos::RCP<CompositeVectorSpace> dep_fac =
          Teuchos::rcp(new CompositeVectorSpace(*my_fac));
      dep_fac->SetOwned(false);

      // Loop over my dependencies, ensuring they meet the requirements.
      for (KeySet::const_iterator key=dependencies_.begin();
           key!=dependencies_.end(); ++key) {
        if (*key != surf_rel_perm_key_) {
          Teuchos::RCP<CompositeVectorSpace> fac = S->RequireField(*key);
          fac->Update(*dep_fac);
        }
      }

      // Recurse into the tree to propagate info to leaves.
      for (KeySet::const_iterator key=dependencies_.begin();
           key!=dependencies_.end(); ++key) {
        S->RequireFieldEvaluator(*key)->EnsureCompatibility(S);
      }

      // Check the dependency for surf rel perm
      S->RequireField(surf_rel_perm_key_)
          ->SetMesh(S->GetMesh("surface"))
          ->AddComponent("cell",AmanziMesh::CELL,1);
      S->RequireFieldEvaluator(surf_rel_perm_key_)->EnsureCompatibility(S);
      
    }
  }
}


void RelPermEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result) {

  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result->Mesh(), -1);
    wrms_->first->Verify();
  }

  // Evaluate k_rel.
  // -- Evaluate the model to calculate krel on cells.
  const Epetra_MultiVector& sat_c = *S->GetFieldData(sat_key_)
      ->ViewComponent("cell",false);
  Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

  int ncells = res_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    int index = (*wrms_->first)[c];
    res_c[0][c] = std::max(wrms_->second[index]->k_relative(sat_c[0][c]), min_val_);
  }

  // -- Potentially evaluate the model on boundary faces as well.
  if (result->HasComponent("boundary_face")) {
    const Epetra_MultiVector& sat_bf = *S->GetFieldData(sat_key_)
        ->ViewComponent("boundary_face",false);
    Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
    const Epetra_Map& face_map = mesh->face_map(false);
  
    // Evaluate the model to calculate krel.
    AmanziMesh::Entity_ID_List cells;
    int nbfaces = res_bf.MyLength();
    for (unsigned int bf=0; bf!=nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      mesh->face_get_cells(f, AmanziMesh::USED, &cells);
      ASSERT(cells.size() == 1);

      int index = (*wrms_->first)[cells[0]];
      res_bf[0][bf] = std::max(wrms_->second[index]->k_relative(sat_bf[0][bf]), min_val_);
    }
  }

  // Patch k_rel with surface rel perm values
  if (is_surf_) {
    const Epetra_MultiVector& surf_kr = *S->GetFieldData(surf_rel_perm_key_)
        ->ViewComponent("cell",false);
    Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = S->GetMesh(surf_mesh_key_);
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
    const Epetra_Map& face_map = mesh->face_map(false);
    
    unsigned int nsurf_cells = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
    for (unsigned int sc=0; sc!=nsurf_cells; ++sc) {
      // need to map from surface quantity on cells to subsurface boundary_face quantity
      AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
      AmanziMesh::Entity_ID bf = vandelay_map.LID(face_map.GID(f));

      res_bf[0][bf] = std::max(surf_kr[0][sc], min_val_);
    }
  }

  // Potentially scale quantities by dens / visc
  if (is_dens_visc_) {
    // -- Scale cells.
    const Epetra_MultiVector& dens_c = *S->GetFieldData(dens_key_)
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& visc_c = *S->GetFieldData(visc_key_)
        ->ViewComponent("cell",false);

    for (unsigned int c=0; c!=ncells; ++c) {
      res_c[0][c] *= dens_c[0][c] / visc_c[0][c];
    }
        
    // Potentially scale boundary faces.
    if (result->HasComponent("boundary_face")) {
      const Epetra_MultiVector& dens_bf = *S->GetFieldData(dens_key_)
          ->ViewComponent("boundary_face",false);
      const Epetra_MultiVector& visc_bf = *S->GetFieldData(visc_key_)
          ->ViewComponent("boundary_face",false);
      Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

      // Evaluate the evaluator to calculate sat.
      int nbfaces = res_bf.MyLength();
      for (unsigned int bf=0; bf!=nbfaces; ++bf) {
        res_bf[0][bf] *= dens_bf[0][bf] / visc_bf[0][bf];
      }
    }
  }

  // Finally, scale by a permeability rescaling from absolute perm.
  result->Scale(1./perm_scale_);
}


void RelPermEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
        Key wrt_key, const Teuchos::Ptr<CompositeVector>& result) {

  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result->Mesh(), -1);
    wrms_->first->Verify();
  }

  if (wrt_key == sat_key_) {
    // dkr / dsl = rho/mu * dkr/dpc * dpc/dsl
    
    // Evaluate k_rel.
    // -- Evaluate the model to calculate krel on cells.
    const Epetra_MultiVector& sat_c = *S->GetFieldData(sat_key_)
        ->ViewComponent("cell",false);
    Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

    int ncells = res_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      int index = (*wrms_->first)[c];
      res_c[0][c] = wrms_->second[index]->d_k_relative(sat_c[0][c]);
      ASSERT(res_c[0][c] >= 0.);
    }

    // -- Potentially evaluate the model on boundary faces as well.
    if (result->HasComponent("boundary_face")) {
      const Epetra_MultiVector& sat_bf = *S->GetFieldData(sat_key_)
          ->ViewComponent("boundary_face",false);
      Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

      Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
      const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
      const Epetra_Map& face_map = mesh->face_map(false);
  
      // Evaluate the model to calculate krel.
      AmanziMesh::Entity_ID_List cells;
      int nbfaces = res_bf.MyLength();
      for (unsigned int bf=0; bf!=nbfaces; ++bf) {
        // given a boundary face, we need the internal cell to choose the right WRM
        AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
        mesh->face_get_cells(f, AmanziMesh::USED, &cells);
        ASSERT(cells.size() == 1);

        int index = (*wrms_->first)[cells[0]];
        res_bf[0][bf] = wrms_->second[index]->d_k_relative(sat_bf[0][bf]);
        ASSERT(res_bf[0][bf] >= 0.);
      }
    }

    // Patch k_rel with surface rel perm values
    if (is_surf_) {
      const Epetra_MultiVector& surf_kr = *S->GetFieldData(surf_rel_perm_key_)
          ->ViewComponent("cell",false);
      Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

      Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = S->GetMesh(surf_mesh_key_);
      Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
      const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
      const Epetra_Map& face_map = mesh->face_map(false);
    
      unsigned int nsurf_cells = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
      for (unsigned int sc=0; sc!=nsurf_cells; ++sc) {
        // need to map from surface quantity on cells to subsurface boundary_face quantity
        AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
        AmanziMesh::Entity_ID bf = vandelay_map.LID(face_map.GID(f));

        //        res_bf[0][bf] = std::max(surf_kr[0][sc], min_val_);
        res_bf[0][bf] = 0.;
      }
    }

    // Potentially scale quantities by dens / visc
    if (is_dens_visc_) {
      // -- Scale cells.
      const Epetra_MultiVector& dens_c = *S->GetFieldData(dens_key_)
          ->ViewComponent("cell",false);
      const Epetra_MultiVector& visc_c = *S->GetFieldData(visc_key_)
          ->ViewComponent("cell",false);

      for (unsigned int c=0; c!=ncells; ++c) {
        res_c[0][c] *= dens_c[0][c] / visc_c[0][c];
      }
        
      // Potentially scale boundary faces.
      if (result->HasComponent("boundary_face")) {
        const Epetra_MultiVector& dens_bf = *S->GetFieldData(dens_key_)
            ->ViewComponent("boundary_face",false);
        const Epetra_MultiVector& visc_bf = *S->GetFieldData(visc_key_)
            ->ViewComponent("boundary_face",false);
        Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

        // Evaluate the evaluator to calculate sat.
        int nbfaces = res_bf.MyLength();
        for (unsigned int bf=0; bf!=nbfaces; ++bf) {
          res_bf[0][bf] *= dens_bf[0][bf] / visc_bf[0][bf];
        }
      }
    }

    // rescale as neeeded
    result->Scale(1./perm_scale_);
    

  } else if (wrt_key == dens_key_) {
    ASSERT(is_dens_visc_);
    // note density > 0
    const Epetra_MultiVector& dens_c = *S->GetFieldData(dens_key_)
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& kr_c = *S->GetFieldData(my_key_)
        ->ViewComponent("cell",false);
    Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

    int ncells = res_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      res_c[0][c] = kr_c[0][c] / dens_c[0][c];
    }

    // Potentially scale boundary faces.
    if (result->HasComponent("boundary_face")) {
      const Epetra_MultiVector& dens_bf = *S->GetFieldData(dens_key_)
          ->ViewComponent("boundary_face",false);
      const Epetra_MultiVector& kr_bf = *S->GetFieldData(my_key_)
          ->ViewComponent("boundary_face",false);
      Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

      // Evaluate the evaluator to calculate sat.
      int nbfaces = res_bf.MyLength();
      for (unsigned int bf=0; bf!=nbfaces; ++bf) {
        res_bf[0][bf] = kr_bf[0][bf] / dens_bf[0][bf];
      }
    }

    
  } else if (wrt_key == visc_key_) {
    ASSERT(is_dens_visc_);
    // note density > 0
    const Epetra_MultiVector& visc_c = *S->GetFieldData(visc_key_)
        ->ViewComponent("cell",false);
    const Epetra_MultiVector& kr_c = *S->GetFieldData(my_key_)
        ->ViewComponent("cell",false);
    Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

    int ncells = res_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      res_c[0][c] = -kr_c[0][c] / visc_c[0][c];
    }


    // Potentially scale boundary faces.
    if (result->HasComponent("boundary_face")) {
      const Epetra_MultiVector& visc_bf = *S->GetFieldData(visc_key_)
          ->ViewComponent("boundary_face",false);
      const Epetra_MultiVector& kr_bf = *S->GetFieldData(my_key_)
          ->ViewComponent("boundary_face",false);
      Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

      // Evaluate the evaluator to calculate sat.
      int nbfaces = res_bf.MyLength();
      for (unsigned int bf=0; bf!=nbfaces; ++bf) {
        res_bf[0][bf] = -kr_bf[0][bf] / visc_bf[0][bf];
      }
    }

  } else {
    ASSERT(0);
  }
  
}



} //namespace
} //namespace
} //namespace
