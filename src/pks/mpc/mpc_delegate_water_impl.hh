/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

// Delegate for heuristic corrections based upon coupled surface/subsurface water.

#pragma once

#include "mpc_surface_subsurface_helpers.hh"

namespace Amanzi {


// Approach 1: global face limiter on the correction size
template <AmanziMesh::Entity_kind FaceEntity>
int
MPCDelegateWater::ModifyCorrection_WaterFaceLimiter(double h,
                                                    Teuchos::RCP<const TreeVector> res,
                                                    Teuchos::RCP<const TreeVector> u,
                                                    Teuchos::RCP<TreeVector> Pu)
{
  int n_modified = 0;
  if (face_limiter_ > 0.) {
    std::string face_entity = AmanziMesh::entity_kind_string(FaceEntity);
    Epetra_MultiVector& domain_Pu_f =
      *Pu->SubVector(i_domain_)->Data()->ViewComponent(face_entity, false);

    for (int f = 0; f != domain_Pu_f.MyLength(); ++f) {
      if (std::abs(domain_Pu_f[0][f]) > face_limiter_) {
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          *vo_->os() << "  LIMITING: dp_old = " << domain_Pu_f[0][f];
        domain_Pu_f[0][f] = domain_Pu_f[0][f] > 0. ? face_limiter_ : -face_limiter_;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          *vo_->os() << ", dp_new = " << domain_Pu_f[0][f] << std::endl;
        n_modified++;
      }
    }
  }
  return n_modified;
}

// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
template <AmanziMesh::Entity_kind FaceEntity>
double
MPCDelegateWater::ModifyCorrection_WaterSpurtDamp(double h,
                                                  Teuchos::RCP<const TreeVector> res,
                                                  Teuchos::RCP<const TreeVector> u,
                                                  Teuchos::RCP<TreeVector> Pu)
{
  // Approach 2
  double damp = 1.;
  if (damp_the_spurt_) {
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    std::string face_entity = AmanziMesh::entity_kind_string(FaceEntity);
    Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(i_domain_)->Data();
    DomainFaceGetter domain_u_f(*domain_mesh_, *domain_u->ViewComponent(face_entity, false));

    Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
    Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);
    DomainFaceGetter domain_Pu_f(*domain_mesh_, *domain_Pu->ViewComponent(face_entity, false));

    int ncells_surf = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int cs = 0; cs != ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
      double p_old = domain_u_f.get<FaceEntity>(f);
      double p_Pu = domain_Pu_f.get<FaceEntity>(f);
      double p_new = p_old - p_Pu;
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        double my_damp = ((patm + cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SPURT (sc=" << surf_mesh_->cell_map(false).GID(cs)
                    << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp
                    << std::endl;
      }
    }

    double proc_damp = damp;
    domain_mesh_->get_comm()->MinAll(&proc_damp, &damp, 1);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SPURT!, coef = " << damp << std::endl;
      domain_Pu->Scale(damp);
    }
  }
  return damp;
}


// Approach 3: capping of the spurt -- limit the max oversaturated pressure
//  if coming from undersaturated.
template <AmanziMesh::Entity_kind FaceEntity>
int
MPCDelegateWater::ModifyCorrection_WaterSpurtCap(double h,
                                                 Teuchos::RCP<const TreeVector> res,
                                                 Teuchos::RCP<const TreeVector> u,
                                                 Teuchos::RCP<TreeVector> Pu,
                                                 double damp)
{
  // Approach 3
  int n_modified = 0;
  if (cap_the_spurt_) {
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    std::string face_entity = AmanziMesh::entity_kind_string(FaceEntity);
    Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(i_domain_)->Data();
    DomainFaceGetter domain_u_f(*domain_mesh_, *domain_u->ViewComponent(face_entity, false));

    Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
    DomainFaceSetter domain_Pu_f(*domain_mesh_, *domain_Pu->ViewComponent(face_entity, false));

    int ncells_surf = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int cs = 0; cs != ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);

      double p_old = domain_u_f.get<FaceEntity>(f);
      double p_Pu = domain_Pu_f.get<FaceEntity>(f);
      double p_new = p_old - p_Pu / damp;
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        double p_corrected = p_old - (patm + cap_size_);
        domain_Pu_f.set<FaceEntity>(f, p_corrected);

        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  CAPPING THE SPURT (sc=" << surf_mesh_->cell_map(false).GID(cs)
                    << ",f=" << u->SubVector(i_domain_)->Data()->Mesh()->face_map(false).GID(f)
                    << "): p_old = " << p_old << ", p_new = " << p_new
                    << ", p_capped = " << p_old - p_corrected << std::endl;
      } else if ((p_new < patm) && (p_old > patm)) {
        // strange attempt to kick NKA when it goes back under?
        // double p_corrected = p_old - (patm - cap_size_);
        // domain_Pu_f.set<FaceEntity>(f, p_corrected);
        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  INVERSE SPURT (sc=" << surf_mesh_->cell_map(false).GID(cs)
                    << "): p_old = " << p_old << ", p_new = " << p_new << std::endl;
      }
    }
  }
  return n_modified;
}


// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
template <AmanziMesh::Entity_kind FaceEntity>
double
MPCDelegateWater::ModifyCorrection_SaturatedSpurtDamp(double h,
                                                      Teuchos::RCP<const TreeVector> res,
                                                      Teuchos::RCP<const TreeVector> u,
                                                      Teuchos::RCP<TreeVector> Pu)
{
  // Approach 2
  double damp = 1.;
  if (damp_the_sat_spurt_) {
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    const Epetra_MultiVector& domain_u_c =
      *u->SubVector(i_domain_)->Data()->ViewComponent("cell", false);
    Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
    Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);

    int ncells_domain =
      domain_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c != ncells_domain; ++c) {
      double p_old = domain_u_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c];
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        double my_damp = ((patm + cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SATURATED SPURT (c=" << domain_mesh_->cell_map(false).GID(c)
                    << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp
                    << std::endl;
      }
    }

    double proc_damp = damp;
    domain_Pu_c.Comm().MinAll(&proc_damp, &damp, 1);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SATURATED SPURT!, coef = " << damp << std::endl;
      domain_Pu->Scale(damp);
    }
  }
  return damp;
}


// Approach 3: capping of the spurt -- limit the max oversaturated pressure
//  if coming from undersaturated.
template <AmanziMesh::Entity_kind FaceEntity>
int
MPCDelegateWater::ModifyCorrection_SaturatedSpurtCap(double h,
                                                     Teuchos::RCP<const TreeVector> res,
                                                     Teuchos::RCP<const TreeVector> u,
                                                     Teuchos::RCP<TreeVector> Pu,
                                                     double damp)
{
  // Approach 3
  int n_modified = 0;
  if (cap_the_sat_spurt_) {
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    const Epetra_MultiVector& domain_u_c =
      *u->SubVector(i_domain_)->Data()->ViewComponent("cell", false);
    Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
    Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);

    int ncells_domain =
      domain_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c != ncells_domain; ++c) {
      double p_old = domain_u_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c] / damp;
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        domain_Pu_c[0][c] = p_old - (patm + cap_size_);
        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  CAPPING THE SATURATED SPURT (c=" << domain_mesh_->cell_map(false).GID(c)
                    << "): p_old = " << p_old << ", p_new = " << p_new
                    << ", p_capped = " << p_old - domain_Pu_c[0][c] << std::endl;
      }
    }
  }

  return n_modified;
}

// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
template <AmanziMesh::Entity_kind FaceEntity>
double
MPCDelegateWater::ModifyCorrection_DesaturatedSpurtDamp(double h,
                                                        Teuchos::RCP<const TreeVector> res,
                                                        Teuchos::RCP<const TreeVector> u,
                                                        Teuchos::RCP<TreeVector> Pu)
{
  // Approach 2
  double damp = 1.;
  if (damp_the_desat_spurt_) {
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    const Epetra_MultiVector& domain_u_c =
      *u->SubVector(i_domain_)->Data()->ViewComponent("cell", false);
    Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
    Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);

    int ncells_domain =
      domain_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c != ncells_domain; ++c) {
      double p_old = domain_u_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c];
      if ((p_new < patm - cap_size_) && (p_old >= patm)) {
        double my_damp = ((patm - cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE DESATURATED SPURT (c="
                    << domain_mesh_->cell_map(false).GID(c) << "): p_old = " << p_old
                    << ", p_new = " << p_new << ", coef = " << my_damp << std::endl;
      }
    }

    double proc_damp = damp;
    domain_Pu_c.Comm().MinAll(&proc_damp, &damp, 1);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE DESATURATED SPURT!, coef = " << damp << std::endl;
      domain_Pu->Scale(damp);
    }
  }
  return damp;
}


// Approach 3: capping of the spurt -- limit the max oversaturated pressure
//  if coming from undersaturated.
template <AmanziMesh::Entity_kind FaceEntity>
int
MPCDelegateWater::ModifyCorrection_DesaturatedSpurtCap(double h,
                                                       Teuchos::RCP<const TreeVector> res,
                                                       Teuchos::RCP<const TreeVector> u,
                                                       Teuchos::RCP<TreeVector> Pu,
                                                       double damp)
{
  // Approach 3
  int n_modified = 0;
  if (cap_the_desat_spurt_) {
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    const Epetra_MultiVector& domain_u_c =
      *u->SubVector(i_domain_)->Data()->ViewComponent("cell", false);
    Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
    Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);

    int ncells_domain =
      domain_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c != ncells_domain; ++c) {
      double p_old = domain_u_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c] / damp;
      if ((p_new < patm - cap_size_) && (p_old >= patm)) {
        domain_Pu_c[0][c] = p_old - (patm - cap_size_);
        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  CAPPING THE DESATURATED SPURT (c=" << domain_mesh_->cell_map(false).GID(c)
                    << "): p_old = " << p_old << ", p_new = " << p_new
                    << ", p_capped = " << p_old - domain_Pu_c[0][c] << std::endl;
      }
    }
  }
  return n_modified;
}


// modify predictor via heuristic stops spurting in the surface flow
template <AmanziMesh::Entity_kind FaceEntity>
bool
MPCDelegateWater::ModifyPredictor_Heuristic(double h, const Teuchos::RCP<TreeVector>& u)
{
  bool modified = false;
  if (modify_predictor_heuristic_) {
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "  MPCWaterCoupler: Modifying predictor with water heuristic" << std::endl;

    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    Teuchos::RCP<CompositeVector> domain_u = u->SubVector(i_domain_)->Data();
    DomainFaceSetter domain_u_f(
      *domain_mesh_, *domain_u->ViewComponent(AmanziMesh::entity_kind_string(FaceEntity), false));

    Epetra_MultiVector& surf_u_c = *u->SubVector(i_surf_)->Data()->ViewComponent("cell", false);

    Key p_surf_key = Keys::getKey(domain_surf_, "pressure");
    const Epetra_MultiVector& surf_u_prev_c =
      *S_->Get<CompositeVector>(p_surf_key, tag_current_).ViewComponent("cell", false);

    int ncells_surf = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int c = 0; c != ncells_surf; ++c) {
      double dp = surf_u_c[0][c] - surf_u_prev_c[0][c];
      double pnew = surf_u_c[0][c] - patm;
      double pold = surf_u_prev_c[0][c] - patm;

      if (pnew > 0) {
        if (dp > pnew) {
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            *vo_->os() << "CHANGING (first over?): p = " << surf_u_c[0][c] << " to "
                       << patm + cap_size_ << std::endl;
          surf_u_c[0][c] = patm + cap_size_;

          AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
          domain_u_f.set<FaceEntity>(f, surf_u_c[0][c]);

        } else if (pold > 0 && dp > pold) {
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            *vo_->os() << "CHANGING (second over?): p = " << surf_u_c[0][c] << " to "
                       << patm + 2 * pold << std::endl;
          surf_u_c[0][c] = patm + 2 * pold;

          AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
          domain_u_f.set<FaceEntity>(f, surf_u_c[0][c]);
        }
      }
    }
    modified = true;
  }

  return modified;
}

// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
// Approach 3: capping of the spurt -- limit the max oversaturated pressure
//  if coming from undersaturated.
template <AmanziMesh::Entity_kind FaceEntity>
bool
MPCDelegateWater::ModifyPredictor_WaterSpurtDamp(double h, const Teuchos::RCP<TreeVector>& u)
{
  // Approach 2
  if (modify_predictor_spurt_damping_) {
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    Epetra_MultiVector& surf_pnew_c = *u->SubVector(i_surf_)->Data()->ViewComponent("cell", false);
    Teuchos::RCP<CompositeVector> domain_pnew = u->SubVector(i_domain_)->Data();
    DomainFaceSetter domain_pnew_f(
      *domain_mesh_,
      *domain_pnew->ViewComponent(AmanziMesh::entity_kind_string(FaceEntity), false));

    Key key_ss = Keys::getKey(domain_ss_, "pressure");
    Teuchos::RCP<const CompositeVector> domain_pold =
      S_->GetPtr<CompositeVector>(key_ss, tag_current_);
    DomainFaceGetter domain_pold_f(
      *domain_mesh_,
      *domain_pold->ViewComponent(AmanziMesh::entity_kind_string(FaceEntity), false));

    double damp = 1.;
    int ncells_surf = surf_mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
    for (int cs = 0; cs != ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
      double p_old = domain_pold_f.get<FaceEntity>(f);
      double p_new = domain_pnew_f.get<FaceEntity>(f);
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        // first over
        double my_damp = ((patm + cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SPURT (1st over) (sc=" << surf_mesh_->cell_map(false).GID(cs)
                    << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp
                    << std::endl;
      } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
        // second over
        double my_damp = ((patm + 2 * (p_old - patm)) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SPURT (2nd over) (sc=" << surf_mesh_->cell_map(false).GID(cs)
                    << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp
                    << std::endl;
      }
    }

    double proc_damp = damp;
    domain_mesh_->get_comm()->MinAll(&proc_damp, &damp, 1);

    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SPURT!, coef = " << damp << std::endl;

      // apply the damping
      db_->WriteVector("p_old", domain_pold.ptr());
      db_->WriteVector("p_new", domain_pnew.ptr());
      domain_pnew->Update(1. - damp, *domain_pold, damp);
      db_->WriteVector("p_damped", domain_pnew.ptr());

      // undamp and cap the surface
      for (unsigned int cs = 0; cs != ncells_surf; ++cs) {
        AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, cs);
        double p_old = domain_pold_f.get<FaceEntity>(f); // THIS IS THE BUG <<<------
        double p_new = domain_pnew_f.get<FaceEntity>(f);

        p_new = (p_new - p_old) / damp + p_old;
        if ((p_new > patm + cap_size_) && (p_old < patm)) {
          // first over
          double new_value = patm + cap_size_;
          domain_pnew_f.set<FaceEntity>(f, new_value);
          surf_pnew_c[0][cs] = new_value;
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            std::cout << "  CAPPING THE SPURT (1st over) (sc="
                      << surf_mesh_->cell_map(false).GID(cs) << "): p_old = " << p_old
                      << ", p_new = " << p_new << ", p_capped = " << new_value << std::endl;
        } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
          // second over
          double new_value = patm + 2 * (p_old - patm);
          domain_pnew_f.set<FaceEntity>(f, new_value);
          surf_pnew_c[0][cs] = new_value;

          if (vo_->os_OK(Teuchos::VERB_HIGH))
            std::cout << "  CAPPING THE SPURT (2nd over) (sc="
                      << surf_mesh_->cell_map(false).GID(cs) << "): p_old = " << p_old
                      << ", p_new = " << p_new << ", p_capped = " << new_value << std::endl;
        } else {
          surf_pnew_c[0][cs] = domain_pnew_f.get<FaceEntity>(f);
        }
      }
    }

    return damp < 1.0;
  }

  return false;
}


// modify predictor via heuristic stops spurting in the surface flow
template <AmanziMesh::Entity_kind FaceEntity>
bool
MPCDelegateWater::ModifyPredictor_TempFromSource(double h, const Teuchos::RCP<TreeVector>& u)
{
  bool modified = false;
  if (modify_predictor_tempfromsource_) {
    modified = true;
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os()
        << "  MPCWaterCoupler: Modifying predictor, taking surface temperature from source."
        << std::endl;

    Epetra_MultiVector& surf_Tnew_c = *u->SubVector(i_Tsurf_)->Data()->ViewComponent("cell", false);
    Epetra_MultiVector& surf_pnew_c = *u->SubVector(i_surf_)->Data()->ViewComponent("cell", false);


    // Epetra_MultiVector& domain_pnew_f = *u->SubVector(i_domain_)->Data()
    //     ->ViewComponent("face",false);
    Teuchos::RCP<CompositeVector> domain_pnew = u->SubVector(i_domain_)->Data();
    DomainFaceSetter domain_pnew_f(
      *domain_mesh_,
      *domain_pnew->ViewComponent(AmanziMesh::entity_kind_string(FaceEntity), false));

    const Epetra_MultiVector& Told =
      *S_->GetPtr<CompositeVector>("surface-temperature", tag_current_)
         ->ViewComponent("cell", false);
    const Epetra_MultiVector& Tsource =
      *S_->GetPtr<CompositeVector>("surface-water_source_temperature", tag_next_)
         ->ViewComponent("cell", false);
    const Epetra_MultiVector& hold =
      *S_->GetPtr<CompositeVector>("surface-ponded_depth", tag_current_)
         ->ViewComponent("cell", false);
    const Epetra_MultiVector& dhsource =
      *S_->GetPtr<CompositeVector>("surface-water_source", tag_current_)
         ->ViewComponent("cell", false);

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->SubVector(i_surf_)->Data()->Mesh();

    for (unsigned int c = 0; c != surf_Tnew_c.MyLength(); ++c) {
      if (surf_Tnew_c[0][c] < 271.15) {
        // frozen, modify predictor to ensure surface is ready to accept ice
        if (surf_pnew_c[0][c] < 101325.) {
          surf_pnew_c[0][c] = 101325.1;
          AmanziMesh::Entity_ID f = surf_mesh_->entity_get_parent(AmanziMesh::CELL, c);
          domain_pnew_f.set<FaceEntity>(f, surf_pnew_c[0][c]);
        }
      }
    }
  }
  return modified;
}

} // namespace Amanzi
