/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! RelPermSutraIceEvaluator: evaluates relative permeability using water retention models.
#include "rel_perm_sutraice_evaluator.hh"
#include "rel_perm_sutraice_drag_term.hh"

namespace Amanzi {
namespace Flow {

// Constructor from ParameterList
RelPermSutraIceEvaluator::RelPermSutraIceEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotypeCV(plist), min_val_(0.)
{
  AMANZI_ASSERT(plist_.isSublist("WRM parameters"));
  Teuchos::ParameterList sublist = plist_.sublist("WRM parameters");
  wrms_ = createWRMPartition(sublist);
  InitializeFromPlist_();
}


RelPermSutraIceEvaluator::RelPermSutraIceEvaluator(Teuchos::ParameterList& plist,
                                                   const Teuchos::RCP<WRMPartition>& wrms)
  : EvaluatorSecondaryMonotypeCV(plist), wrms_(wrms), min_val_(0.)
{
  InitializeFromPlist_();
}

Teuchos::RCP<Evaluator>
RelPermSutraIceEvaluator::Clone() const
{
  return Teuchos::rcp(new RelPermSutraIceEvaluator(*this));
}


void
RelPermSutraIceEvaluator::InitializeFromPlist_()
{
  // dependencies
  Key domain_name = Keys::getDomain(my_keys_.front().first);
  Tag tag = my_keys_.front().second;

  // -- saturation liquid
  sat_key_ = Keys::readKey(plist_, domain_name, "saturation", "saturation_liquid");
  dependencies_.insert(KeyTag{ sat_key_, tag });
  sat_gas_key_ = Keys::readKey(plist_, domain_name, "saturation gas", "saturation_gas");
  dependencies_.insert(KeyTag{ sat_gas_key_, tag });

  is_dens_visc_ = plist_.get<bool>("use density on viscosity in rel perm", true);
  if (is_dens_visc_) {
    dens_key_ = Keys::readKey(plist_, domain_name, "density", "molar_density_liquid");
    dependencies_.insert(KeyTag{ dens_key_, tag });

    visc_key_ = Keys::readKey(plist_, domain_name, "viscosity", "viscosity_liquid");
    dependencies_.insert(KeyTag{ visc_key_, tag });
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
    boundary_krel_ = BoundaryRelPerm::BOUNDARY_PRESSURE;
  } else if (boundary_krel == "interior pressure") {
    boundary_krel_ = BoundaryRelPerm::INTERIOR_PRESSURE;
  } else if (boundary_krel == "harmonic mean") {
    boundary_krel_ = BoundaryRelPerm::HARMONIC_MEAN;
  } else if (boundary_krel == "arithmetic mean") {
    boundary_krel_ = BoundaryRelPerm::ARITHMETIC_MEAN;
  } else if (boundary_krel == "one") {
    boundary_krel_ = BoundaryRelPerm::ONE;
  } else if (boundary_krel == "surface rel perm") {
    boundary_krel_ = BoundaryRelPerm::SURF_REL_PERM;
  } else {
    Errors::Message msg("RelPermEvaluator: parameter \"boundary rel perm strategy\" not valid: "
                        "valid are \"boundary pressure\", \"interior pressure\", \"harmonic "
                        "mean\", \"arithmetic mean\", \"one\", and \"surface rel perm\"");
    throw(msg);
  }

  // surface alterations
  if (boundary_krel_ == BoundaryRelPerm::SURF_REL_PERM) {
    surf_domain_ = Keys::readDomainHint(plist_, domain_name, "domain", "surface");
    surf_rel_perm_key_ = Keys::readKey(plist_,
                                       surf_domain_,
                                       "surface relative permeability",
                                       Keys::getVarName(my_keys_.front().first));
    dependencies_.insert(KeyTag{ surf_rel_perm_key_, tag });
  }

  // cutoff above 0?
  min_val_ = plist_.get<double>("minimum rel perm cutoff", 0.);
  perm_scale_ = plist_.get<double>("permeability rescaling");
  omega_ = plist_.get<double>("coefficient in drag term of kr", 1);
}


// Special purpose EnsureCompatibility required because of surface rel perm.
void
RelPermSutraIceEvaluator::EnsureCompatibility_ToDeps_(State& S)
{
  Key my_key = my_keys_.front().first;
  Tag tag = my_keys_.front().second;
  const auto& my_fac = S.Require<CompositeVector, CompositeVectorSpace>(my_key, tag);

  // If my requirements have not yet been set, we'll have to hope they
  // get set by someone later.  For now just defer.
  if (my_fac.Mesh() != Teuchos::null) {
    // Loop over my dependencies, ensuring they meet the requirements.
    for (const auto& dep : dependencies_) {
      if (dep.first != surf_rel_perm_key_) {
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(my_fac);
      } else {
        CompositeVectorSpace surf_kr_fac;
        surf_kr_fac.SetMesh(S.GetMesh(Keys::getDomain(dep.first)))
          ->AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
        S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(surf_kr_fac);
      }
    }
  }
}


void
RelPermSutraIceEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& result)
{
  result[0]->PutScalar(0.);

  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result[0]->Mesh(), -1);
    wrms_->first->Verify();
  }

  // Evaluate k_rel.
  // -- Evaluate the model to calculate krel on cells.
  Tag tag = my_keys_.front().second;
  const Epetra_MultiVector& sat_c =
    *S.GetPtr<CompositeVector>(sat_key_, tag)->ViewComponent("cell", false);
  const Epetra_MultiVector& sat_gas_c =
    *S.GetPtr<CompositeVector>(sat_gas_key_, tag)->ViewComponent("cell", false);
  Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

  int ncells = res_c.MyLength();
  for (unsigned int c = 0; c != ncells; ++c) {
    int index = (*wrms_->first)[c];

    double sat_res = wrms_->second[index]->residualSaturation();
    double coef = SutraIceTerm::dragcoef(sat_c[0][c], sat_gas_c[0][c], sat_res, omega_);
    res_c[0][c] = std::max(wrms_->second[index]->k_relative(1. - sat_gas_c[0][c]) * coef, min_val_);
  }

  // -- Potentially evaluate the model on boundary faces as well.
  if (result[0]->HasComponent("boundary_face")) {
    const Epetra_MultiVector& sat_bf =
      *S.GetPtr<CompositeVector>(sat_key_, tag)->ViewComponent("boundary_face", false);
    const Epetra_MultiVector& sat_gas_bf =
      *S.GetPtr<CompositeVector>(sat_gas_key_, tag)->ViewComponent("boundary_face", false);
    Epetra_MultiVector& res_bf = *result[0]->ViewComponent("boundary_face", false);

    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result[0]->Mesh();
    const Epetra_Map& vandelay_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);
    const Epetra_Map& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);

    // Evaluate the model to calculate krel.
    int nbfaces = res_bf.MyLength();
    for (unsigned int bf = 0; bf != nbfaces; ++bf) {
      // given a boundary face, we need the internal cell to choose the right WRM
      AmanziMesh::Entity_ID f = face_map.LID(vandelay_map.GID(bf));
      auto cells = mesh->getFaceCells(f, AmanziMesh::Parallel_kind::ALL);
      AMANZI_ASSERT(cells.size() == 1);

      int index = (*wrms_->first)[cells[0]];
      double sat_res = wrms_->second[index]->residualSaturation();
      double krel;
      double coef_b = SutraIceTerm::dragcoef(sat_bf[0][bf], sat_gas_bf[0][bf], sat_res, omega_);
      double coef_c =
        SutraIceTerm::dragcoef(sat_c[0][cells[0]], sat_gas_c[0][cells[0]], sat_res, omega_);

      if (boundary_krel_ == BoundaryRelPerm::HARMONIC_MEAN) {
        double krelb =
          std::max(wrms_->second[index]->k_relative(1. - sat_gas_bf[0][bf]) * coef_b, min_val_);
        double kreli = std::max(
          wrms_->second[index]->k_relative(1. - sat_gas_c[0][cells[0]]) * coef_c, min_val_);
        krel = 1.0 / (1.0 / krelb + 1.0 / kreli);
      } else if (boundary_krel_ == BoundaryRelPerm::ARITHMETIC_MEAN) {
        double krelb =
          std::max(wrms_->second[index]->k_relative(1. - sat_gas_bf[0][bf]) * coef_b, min_val_);
        double kreli = std::max(
          wrms_->second[index]->k_relative(1. - sat_gas_c[0][cells[0]]) * coef_c, min_val_);
        krel = (krelb + kreli) / 2.0;
      } else if (boundary_krel_ == BoundaryRelPerm::INTERIOR_PRESSURE) {
        krel = std::max(wrms_->second[index]->k_relative(1. - sat_gas_c[0][cells[0]]) * coef_c,
                        min_val_);
      } else if (boundary_krel_ == BoundaryRelPerm::ONE) {
        krel = 1.;
      } else {
        krel = wrms_->second[index]->k_relative(1. - sat_gas_bf[0][bf]) * coef_b;
      }
      res_bf[0][bf] = std::max(krel, min_val_);
    }
  }

  // Patch k_rel with surface rel perm values
  if (boundary_krel_ == BoundaryRelPerm::SURF_REL_PERM) {
    const Epetra_MultiVector& surf_kr =
      *S.GetPtr<CompositeVector>(surf_rel_perm_key_, tag)->ViewComponent("cell", false);
    Epetra_MultiVector& res_bf = *result[0]->ViewComponent("boundary_face", false);

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = S.GetMesh(surf_domain_);
    Teuchos::RCP<const AmanziMesh::Mesh> mesh = result[0]->Mesh();
    const Epetra_Map& vandelay_map = mesh->getMap(AmanziMesh::Entity_kind::BOUNDARY_FACE,false);
    const Epetra_Map& face_map = mesh->getMap(AmanziMesh::Entity_kind::FACE,false);

    unsigned int nsurf_cells =
      surf_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);
    for (unsigned int sc = 0; sc != nsurf_cells; ++sc) {
      // need to map from surface quantity on cells to subsurface boundary_face quantity
      AmanziMesh::Entity_ID f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, sc);
      AmanziMesh::Entity_ID bf = vandelay_map.LID(face_map.GID(f));

      res_bf[0][bf] = std::max(surf_kr[0][sc], min_val_);
    }
  }

  // Potentially scale quantities by dens / visc
  if (is_dens_visc_) {
    // -- Scale cells.
    const Epetra_MultiVector& dens_c =
      *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("cell", false);
    const Epetra_MultiVector& visc_c =
      *S.GetPtr<CompositeVector>(visc_key_, tag)->ViewComponent("cell", false);

    for (unsigned int c = 0; c != ncells; ++c) { res_c[0][c] *= dens_c[0][c] / visc_c[0][c]; }

    // Potentially scale boundary faces.
    if (result[0]->HasComponent("boundary_face")) {
      const Epetra_MultiVector& dens_bf =
        *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("boundary_face", false);
      const Epetra_MultiVector& visc_bf =
        *S.GetPtr<CompositeVector>(visc_key_, tag)->ViewComponent("boundary_face", false);
      Epetra_MultiVector& res_bf = *result[0]->ViewComponent("boundary_face", false);

      // Evaluate the evaluator to calculate sat.
      int nbfaces = res_bf.MyLength();
      for (unsigned int bf = 0; bf != nbfaces; ++bf) {
        res_bf[0][bf] *= dens_bf[0][bf] / visc_bf[0][bf];
      }
    }
  }

  // Finally, scale by a permeability rescaling from absolute perm.
  result[0]->Scale(1. / perm_scale_);
}


void
RelPermSutraIceEvaluator::EvaluatePartialDerivative_(const State& S,
                                                     const Key& wrt_key,
                                                     const Tag& wrt_tag,
                                                     const std::vector<CompositeVector*>& result)
{
  // Initialize the MeshPartition
  if (!wrms_->first->initialized()) {
    wrms_->first->Initialize(result[0]->Mesh(), -1);
    wrms_->first->Verify();
  }

  Tag tag = my_keys_.front().second;

  if (wrt_key == sat_key_) {
    // dkr / dsl = rho/mu * dkr/dpc * dpc/dsl

    // Evaluate k_rel.
    // -- Evaluate the model to calculate krel on cells.
    const Epetra_MultiVector& sat_c =
      *S.GetPtr<CompositeVector>(sat_key_, tag)->ViewComponent("cell", false);
    const Epetra_MultiVector& sat_gas_c =
      *S.GetPtr<CompositeVector>(sat_gas_key_, tag)->ViewComponent("cell", false);
    Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

    int ncells = res_c.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) {
      int index = (*wrms_->first)[c];
      double sat_res = wrms_->second[index]->residualSaturation();
      double dcoef_dsl =
        SutraIceTerm::d_dragcoef_dsl(sat_c[0][c], sat_gas_c[0][c], sat_res, omega_);
      res_c[0][c] = wrms_->second[index]->k_relative(1. - sat_gas_c[0][c]) * dcoef_dsl;
      AMANZI_ASSERT(res_c[0][c] >= 0.);
    }

    // -- Potentially evaluate the model on boundary faces as well.
    if (result[0]->HasComponent("boundary_face")) {
      // const Epetra_MultiVector& sat_bf = *S.GetPtr<CompositeVector>(sat_key_)
      //                                    ->ViewComponent("boundary_face",false);
      Epetra_MultiVector& res_bf = *result[0]->ViewComponent("boundary_face", false);

      // it is unclear that this is used -- in fact it probably isn't -- so we
      // just set it to zero --etc
      res_bf.PutScalar(0.);
    }

    // Potentially scale quantities by dens / visc
    if (is_dens_visc_) {
      // -- Scale cells.
      const Epetra_MultiVector& dens_c =
        *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("cell", false);
      const Epetra_MultiVector& visc_c =
        *S.GetPtr<CompositeVector>(visc_key_, tag)->ViewComponent("cell", false);

      for (unsigned int c = 0; c != ncells; ++c) { res_c[0][c] *= dens_c[0][c] / visc_c[0][c]; }
    }

    // rescale as neeeded
    result[0]->Scale(1. / perm_scale_);

  } else if (wrt_key == sat_gas_key_) {
    // Evaluate k_rel.
    // -- Evaluate the model to calculate krel on cells.
    const Epetra_MultiVector& sat_c =
      *S.GetPtr<CompositeVector>(sat_key_, tag)->ViewComponent("cell", false);
    const Epetra_MultiVector& sat_gas_c =
      *S.GetPtr<CompositeVector>(sat_gas_key_, tag)->ViewComponent("cell", false);
    Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

    int ncells = res_c.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) {
      int index = (*wrms_->first)[c];
      double sat_res = wrms_->second[index]->residualSaturation();
      double coef = SutraIceTerm::dragcoef(sat_c[0][c], sat_gas_c[0][c], sat_res, omega_);
      double dcoef_dsg =
        SutraIceTerm::d_dragcoef_dsg(sat_c[0][c], sat_gas_c[0][c], sat_res, omega_);
      res_c[0][c] = wrms_->second[index]->d_k_relative(1. - sat_gas_c[0][c]) * (-1) * coef +
                    wrms_->second[index]->k_relative(1. - sat_gas_c[0][c]) * dcoef_dsg;
    }

    // -- Potentially evaluate the model on boundary faces as well.
    if (result[0]->HasComponent("boundary_face")) {
      // const Epetra_MultiVector& sat_bf = *S.GetPtr<CompositeVector>(sat_key_)
      //                                    ->ViewComponent("boundary_face",false);
      Epetra_MultiVector& res_bf = *result[0]->ViewComponent("boundary_face", false);

      // it is unclear that this is used -- in fact it probably isn't -- so we
      // just set it to zero --etc
      res_bf.PutScalar(0.);
    }

    // Potentially scale quantities by dens / visc
    if (is_dens_visc_) {
      // -- Scale cells.
      const Epetra_MultiVector& dens_c =
        *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("cell", false);
      const Epetra_MultiVector& visc_c =
        *S.GetPtr<CompositeVector>(visc_key_, tag)->ViewComponent("cell", false);

      for (unsigned int c = 0; c != ncells; ++c) { res_c[0][c] *= dens_c[0][c] / visc_c[0][c]; }
    }

    // rescale as neeeded
    result[0]->Scale(1. / perm_scale_);

  } else if (wrt_key == dens_key_) {
    AMANZI_ASSERT(is_dens_visc_);
    // note density > 0
    const Epetra_MultiVector& dens_c =
      *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("cell", false);
    const Epetra_MultiVector& kr_c =
      *S.GetPtr<CompositeVector>(my_keys_.front().first, tag)->ViewComponent("cell", false);
    Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

    int ncells = res_c.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) { res_c[0][c] = kr_c[0][c] / dens_c[0][c]; }

    // Potentially scale boundary faces.
    if (result[0]->HasComponent("boundary_face")) {
      const Epetra_MultiVector& dens_bf =
        *S.GetPtr<CompositeVector>(dens_key_, tag)->ViewComponent("boundary_face", false);
      const Epetra_MultiVector& kr_bf = *S.GetPtr<CompositeVector>(my_keys_.front().first, tag)
                                           ->ViewComponent("boundary_face", false);
      Epetra_MultiVector& res_bf = *result[0]->ViewComponent("boundary_face", false);

      // Evaluate the evaluator to calculate sat.
      int nbfaces = res_bf.MyLength();
      for (unsigned int bf = 0; bf != nbfaces; ++bf) {
        res_bf[0][bf] = kr_bf[0][bf] / dens_bf[0][bf];
      }
    }


  } else if (wrt_key == visc_key_) {
    AMANZI_ASSERT(is_dens_visc_);
    // note density > 0
    const Epetra_MultiVector& visc_c =
      *S.GetPtr<CompositeVector>(visc_key_, tag)->ViewComponent("cell", false);
    const Epetra_MultiVector& kr_c =
      *S.GetPtr<CompositeVector>(my_keys_.front().first, tag)->ViewComponent("cell", false);
    Epetra_MultiVector& res_c = *result[0]->ViewComponent("cell", false);

    int ncells = res_c.MyLength();
    for (unsigned int c = 0; c != ncells; ++c) { res_c[0][c] = -kr_c[0][c] / visc_c[0][c]; }


    // Potentially scale boundary faces.
    if (result[0]->HasComponent("boundary_face")) {
      const Epetra_MultiVector& visc_bf =
        *S.GetPtr<CompositeVector>(visc_key_, tag)->ViewComponent("boundary_face", false);
      const Epetra_MultiVector& kr_bf = *S.GetPtr<CompositeVector>(my_keys_.front().first, tag)
                                           ->ViewComponent("boundary_face", false);
      Epetra_MultiVector& res_bf = *result[0]->ViewComponent("boundary_face", false);

      // Evaluate the evaluator to calculate sat.
      int nbfaces = res_bf.MyLength();
      for (unsigned int bf = 0; bf != nbfaces; ++bf) {
        res_bf[0][bf] = -kr_bf[0][bf] / visc_bf[0][bf];
      }
    }

  } else {
    AMANZI_ASSERT(0);
  }
}

} // namespace Flow
} // namespace Amanzi
