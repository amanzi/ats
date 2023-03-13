/*
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/
//! RelPermSutraIceEvaluator: evaluates relative permeability using water retention models.

#include "rel_perm_frzcampbell_evaluator.hh"

namespace Amanzi {
namespace Flow {

// Constructor from ParameterList
RelPermFrzCampbellEvaluator::RelPermFrzCampbellEvaluator(Teuchos::ParameterList& plist) :
    SecondaryVariableFieldEvaluator(plist),
    min_val_(0.) {

  AMANZI_ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(sublist);
  InitializeFromPlist_();
}


RelPermFrzCampbellEvaluator::RelPermFrzCampbellEvaluator(Teuchos::ParameterList& plist,
        const Teuchos::RCP<WRMPartition>& wrms) :
    SecondaryVariableFieldEvaluator(plist),
    wrms_(wrms),
    min_val_(0.) {
  InitializeFromPlist_();
}

Teuchos::RCP<FieldEvaluator>
RelPermFrzCampbellEvaluator::Clone() const {
  return Teuchos::rcp(new RelPermFrzCampbellEvaluator(*this));
}


void RelPermFrzCampbellEvaluator::InitializeFromPlist_() {
  // my keys are for saturation and rel perm.
  if (my_key_ == std::string("")) {
    my_key_ = plist_.get<std::string>("rel perm key", "relative_permeability");
  }

  // dependencies
  Key domain_name = Keys::getDomain(my_key_);

  // -- saturation liquid
  sat_key_ = Keys::readKey(plist_, domain_name, "saturation", "saturation_liquid");
  dependencies_.insert(sat_key_);

  // -- saturation gas
  sat_gas_key_ = Keys::readKey(plist_, domain_name, "saturation gas", "saturation_gas");
  dependencies_.insert(sat_gas_key_);

  is_dens_visc_ = plist_.get<bool>("use density on viscosity in rel perm", true);
  if (is_dens_visc_) {
    dens_key_ = Keys::readKey(plist_, domain_name, "density", "molar_density_liquid");
    dependencies_.insert(dens_key_);

    visc_key_ = Keys::readKey(plist_, domain_name, "viscosity", "viscosity_liquid");
    dependencies_.insert(visc_key_);
  }

  // boundary rel perm settings -- deals with deprecated option
  std::string boundary_krel = "boundary pressure";
  if (plist_.isParameter("boundary rel perm strategy")) {
    boundary_krel = plist_.get<std::string>("boundary rel perm strategy", "boundary pressure");
  } else if (plist_.isParameter("use surface rel perm") &&
             plist_.get<bool>("use surface rel perm")) {
    boundary_krel = "surface rel perm";
  }

  if (boundary_krel == "boundary pressure") {
    boundary_krel_ = BoundaryFrzCampbellRelPerm::BOUNDARY_PRESSURE;
  } else if (boundary_krel == "interior pressure") {
    boundary_krel_ = BoundaryFrzCampbellRelPerm::INTERIOR_PRESSURE;
  } else if (boundary_krel == "harmonic mean") {
    boundary_krel_ = BoundaryFrzCampbellRelPerm::HARMONIC_MEAN;
  } else if (boundary_krel == "arithmetic mean") {
    boundary_krel_ = BoundaryFrzCampbellRelPerm::ARITHMETIC_MEAN;
  } else if (boundary_krel == "one") {
    boundary_krel_ = BoundaryFrzCampbellRelPerm::ONE;
  } else if (boundary_krel == "surface rel perm") {
    boundary_krel_ = BoundaryFrzCampbellRelPerm::SURF_REL_PERM;
  } else {
    Errors::Message msg("RelPermFrzCampbellEvaluator: parameter \"boundary rel perm strategy\" not valid: valid are \"boundary pressure\", \"interior pressure\", \"harmonic mean\", \"arithmetic mean\", \"one\", and \"surface rel perm\"");
    throw(msg);
  }

  // surface alterations
  if (boundary_krel_ == BoundaryFrzCampbellRelPerm::SURF_REL_PERM) {
    surf_domain_ = Keys::readDomainHint(plist_, domain_name, "domain", "surface");
    surf_rel_perm_key_ = Keys::readKey(plist_, surf_domain_, "surface relative permeability", Keys::getVarName(my_key_));
    dependencies_.insert(surf_rel_perm_key_);
  }

  // cutoff above 0?
  min_val_ = plist_.get<double>("minimum rel perm cutoff", 0.);
  perm_scale_ = plist_.get<double>("permeability rescaling");
  omega_ = plist_.get<double>("scale dependent parameter", 3.0);
  b_ = plist_.get<double>("clapp hornberger b", 2.0);
}


// Special purpose EnsureCompatibility required because of surface rel perm.
void RelPermFrzCampbellEvaluator::EnsureCompatibility(const Teuchos::Ptr<State>& S) {
  if (boundary_krel_ != BoundaryFrzCampbellRelPerm::SURF_REL_PERM) {
    SecondaryVariableFieldEvaluator::EnsureCompatibility(S);
  } else {

    // Ensure my field exists.  Requirements should be already set.
    AMANZI_ASSERT(my_key_ != std::string(""));
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

      Key domain = Keys::getDomain(surf_rel_perm_key_);
      S->RequireField(surf_rel_perm_key_)
          ->SetMesh(S->GetMesh(domain))
          ->AddComponent("cell",AmanziMesh::CELL,1);

      S->RequireFieldEvaluator(surf_rel_perm_key_)->EnsureCompatibility(S);

    }
  }
}


void RelPermFrzCampbellEvaluator::EvaluateField_(const Teuchos::Ptr<State>& S,
        const Teuchos::Ptr<CompositeVector>& result)
{
  result->PutScalar(0.);

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

  const Epetra_MultiVector& sat_gas_c = *S->GetFieldData(sat_gas_key_)
      ->ViewComponent("cell",false);


  int ncells = res_c.MyLength();
  for (unsigned int c=0; c!=ncells; ++c) {
    double coef = 1 - std::exp(-omega_*(sat_c[0][c]+sat_gas_c[0][c])) + std::exp(-omega_);
    res_c[0][c] = std::max(coef*std::pow(1-sat_gas_c[0][c], 3+2*b_), min_val_);
  }

  // -- Potentially evaluate the model on boundary faces as well.
  if (result->HasComponent("boundary_face")) {
    const Epetra_MultiVector& sat_bf = *S->GetFieldData(sat_key_)
                                       ->ViewComponent("boundary_face",false);
    Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

    const Epetra_MultiVector& sat_gas_bf = *S->GetFieldData(sat_gas_key_)
                                           ->ViewComponent("boundary_face",false);


    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
    const Epetra_Map& face_map = mesh->face_map(false);

    // Evaluate the model to calculate krel.
    AmanziMesh::Entity_ID_List cells;
    int nbfaces = res_bf.MyLength();
    for (unsigned int bf=0; bf!=nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      mesh->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      double krel;
      double coef_b = 1 - std::exp(-omega_*(sat_bf[0][bf]+sat_gas_bf[0][bf])) + std::exp(-omega_); 
      double coef_c = 1 - std::exp(-omega_*(sat_c[0][cells[0]]+sat_gas_c[0][cells[0]])) + std::exp(-omega_);
      if (boundary_krel_ == BoundaryFrzCampbellRelPerm::HARMONIC_MEAN) {
        double krelb = std::max(coef_b*std::pow(1-sat_gas_bf[0][bf], 3+2*b_), min_val_);
        double kreli = std::max(coef_c*std::pow(1-sat_gas_c[0][cells[0]], 3+2*b_), min_val_);
        krel = 1.0 / (1.0/krelb + 1.0/kreli);
      } else if (boundary_krel_ == BoundaryFrzCampbellRelPerm::ARITHMETIC_MEAN) {
        double krelb = std::max(coef_b*std::pow(1-sat_gas_bf[0][bf], 3+2*b_), min_val_);
        double kreli = std::max(coef_c*std::pow(1-sat_gas_c[0][cells[0]], 3+2*b_), min_val_);
        krel = (krelb + kreli)/2.0;
      } else if (boundary_krel_ == BoundaryFrzCampbellRelPerm::INTERIOR_PRESSURE) {
        krel = std::max(coef_c*std::pow(1-sat_gas_c[0][cells[0]], 3+2*b_), min_val_);
      } else if (boundary_krel_ == BoundaryFrzCampbellRelPerm::ONE) {
        krel = 1.;
      } else {
        krel = coef_b*std::pow(1-sat_gas_bf[0][bf], 3+2*b_);
      }
      res_bf[0][bf] = std::max(krel, min_val_);
    }
  }

  // Patch k_rel with surface rel perm values
  if (boundary_krel_ == BoundaryFrzCampbellRelPerm::SURF_REL_PERM) {
    const Epetra_MultiVector& surf_kr = *S->GetFieldData(surf_rel_perm_key_)
        ->ViewComponent("cell",false);
    Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = S->GetMesh(surf_domain_);
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
    const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
    const Epetra_Map& face_map = mesh->face_map(false);

    unsigned int nsurf_cells = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
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


void RelPermFrzCampbellEvaluator::EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S,
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

    const Epetra_MultiVector& sat_gas_c = *S->GetFieldData(sat_gas_key_)
        ->ViewComponent("cell",false);

    int ncells = res_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      res_c[0][c] = omega_*std::pow((1.-sat_gas_c[0][c]), 2*b_+3)*std::exp(-omega_*(sat_c[0][c]+sat_gas_c[0][c]));
      AMANZI_ASSERT(res_c[0][c] >= 0.);
    }

    // -- Potentially evaluate the model on boundary faces as well.
    if (result->HasComponent("boundary_face")) {
      Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

      // it is unclear that this is used -- in fact it probably isn't --etc
      res_bf.PutScalar(0.);
    }

    // Patch k_rel with surface rel perm values
    if (boundary_krel_ == BoundaryFrzCampbellRelPerm::SURF_REL_PERM) {
      const Epetra_MultiVector& surf_kr = *S->GetFieldData(surf_rel_perm_key_)
          ->ViewComponent("cell",false);
      Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

      Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = S->GetMesh(surf_domain_);
      Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
      const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
      const Epetra_Map& face_map = mesh->face_map(false);

      unsigned int nsurf_cells = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
      for (unsigned int sc=0; sc!=nsurf_cells; ++sc) {
        // need to map from surface quantity on cells to subsurface boundary_face quantity
        AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
        AmanziMesh::Entity_ID bf = vandelay_map.LID(face_map.GID(f));

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

  } else if (wrt_key == sat_gas_key_) {

    // Evaluate k_rel.
    // -- Evaluate the model to calculate krel on cells.
    const Epetra_MultiVector& sat_c = *S->GetFieldData(sat_key_)
        ->ViewComponent("cell",false);
    Epetra_MultiVector& res_c = *result->ViewComponent("cell",false);

    const Epetra_MultiVector& sat_gas_c = *S->GetFieldData(sat_gas_key_)
        ->ViewComponent("cell",false);

    int ncells = res_c.MyLength();
    for (unsigned int c=0; c!=ncells; ++c) {
      double p0 = (2*b_+3)*std::pow((1-sat_gas_c[0][c]), 2*b_+2);
      double p1 = std::exp(-omega_*(sat_c[0][c]+sat_gas_c[0][c]))-std::exp(-omega_)-1.;
      double p2 = omega_*std::exp(-omega_*(sat_c[0][c]+sat_gas_c[0][c]))*std::pow(1.-sat_gas_c[0][c], 2*b_+3);
      res_c[0][c] = p0*p1+p2;
    }

    // -- Potentially evaluate the model on boundary faces as well.
    if (result->HasComponent("boundary_face")) {
      Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

      // it is unclear that this is used -- in fact it probably isn't --etc
      res_bf.PutScalar(0.);
    }

    // Patch k_rel with surface rel perm values
    if (boundary_krel_ == BoundaryFrzCampbellRelPerm::SURF_REL_PERM) {
      const Epetra_MultiVector& surf_kr = *S->GetFieldData(surf_rel_perm_key_)
          ->ViewComponent("cell",false);
      Epetra_MultiVector& res_bf = *result->ViewComponent("boundary_face",false);

      Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = S->GetMesh(surf_domain_);
      Teuchos::RCP<const AmanziMesh::Mesh> mesh = result->Mesh();
      const Epetra_Map& vandelay_map = mesh->exterior_face_map(false);
      const Epetra_Map& face_map = mesh->face_map(false);

      unsigned int nsurf_cells = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
      for (unsigned int sc=0; sc!=nsurf_cells; ++sc) {
        // need to map from surface quantity on cells to subsurface boundary_face quantity
        AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, sc);
        AmanziMesh::Entity_ID bf = vandelay_map.LID(face_map.GID(f));

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
    AMANZI_ASSERT(is_dens_visc_);
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
    AMANZI_ASSERT(is_dens_visc_);
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
    AMANZI_ASSERT(0);
  }

}



} //namespace
} //namespace
