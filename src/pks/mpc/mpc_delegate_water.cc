/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

// Delegate for heuristic corrections based upon coupled surface/subsurface water.

#include "mpc_delegate_water.hh"
#include "mpc_surface_subsurface_helpers.hh"

namespace Amanzi {


MPCDelegateWater::MPCDelegateWater(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                                   const Teuchos::RCP<State>& S,
                                   const std::string& domain)
  : plist_(plist),
    S_(S),
    i_domain_(-1),
    i_surf_(-1),
    i_Tdomain_(-1),
    i_Tsurf_(-1),
    domain_ss_(domain)
{
  // predictor control
  modify_predictor_heuristic_ = plist_->get<bool>("modify predictor with heuristic", false);
  modify_predictor_spurt_damping_ =
    plist_->get<bool>("modify predictor damp and cap the water spurt", false);
  modify_predictor_tempfromsource_ =
    plist_->get<bool>("modify predictor surface temperature from source", false);

  // precon control
  face_limiter_ = plist_->get<double>("global water face limiter", -1.0);
  cap_the_spurt_ = plist_->get<bool>("cap the water spurt", false);
  damp_the_spurt_ = plist_->get<bool>("damp the water spurt", false);
  bool damp_and_cap_the_spurt = plist_->get<bool>("damp and cap the water spurt", false);
  if (damp_and_cap_the_spurt) {
    damp_the_spurt_ = true;
    cap_the_spurt_ = true;
  }

  cap_the_sat_spurt_ = plist_->get<bool>("cap the saturated spurt", false);
  damp_the_sat_spurt_ = plist_->get<bool>("damp the saturated spurt", false);
  bool damp_and_cap_the_sat_spurt = plist_->get<bool>("damp and cap the saturated spurt", false);
  if (damp_and_cap_the_sat_spurt) {
    damp_the_sat_spurt_ = true;
    cap_the_sat_spurt_ = true;
  }

  cap_the_desat_spurt_ = plist_->get<bool>("cap the desaturated spurt", false);
  damp_the_desat_spurt_ = plist_->get<bool>("damp the desaturated spurt", false);
  bool damp_and_cap_the_desat_spurt =
    plist_->get<bool>("damp and cap the desaturated spurt", false);
  if (damp_and_cap_the_desat_spurt) {
    damp_the_desat_spurt_ = true;
    cap_the_desat_spurt_ = true;
  }

  // set the size of the caps
  if (cap_the_spurt_ || damp_the_spurt_ || cap_the_sat_spurt_ || damp_the_sat_spurt_ ||
      modify_predictor_heuristic_ || modify_predictor_spurt_damping_) {
    cap_size_ = plist_->get<double>("cap over atmospheric", 100.0);
  }

  // create the VO
  vo_ = Teuchos::rcp(new VerboseObject(plist->name(), *plist_));
}


// Approach 1: global face limiter on the correction size
int
MPCDelegateWater::ModifyCorrection_WaterFaceLimiter(double h,
                                                    Teuchos::RCP<const TreeVector> res,
                                                    Teuchos::RCP<const TreeVector> u,
                                                    Teuchos::RCP<TreeVector> Pu)
{
  int n_modified = 0;
  if (face_limiter_ > 0.) {
    std::string face_entity;
    if (Pu->SubVector(i_domain_)->Data()->HasComponent("face")) {
      face_entity = "face";
    } else if (Pu->SubVector(i_domain_)->Data()->HasComponent("boundary_face")) {
      face_entity = "boundary_face";
    } else {
      Errors::Message message("Subsurface vector does not have face component.");
      Exceptions::amanzi_throw(message);
    }

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
double
MPCDelegateWater::ModifyCorrection_WaterSpurtDamp(double h,
                                                  Teuchos::RCP<const TreeVector> res,
                                                  Teuchos::RCP<const TreeVector> u,
                                                  Teuchos::RCP<TreeVector> Pu)
{
  const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  std::string face_entity;
  if (Pu->SubVector(i_domain_)->Data()->HasComponent("face")) {
    face_entity = "face";
  } else if (Pu->SubVector(i_domain_)->Data()->HasComponent("boundary_face")) {
    face_entity = "boundary_face";
  } else {
    Errors::Message message("Subsurface vector does not have face component.");
    Exceptions::amanzi_throw(message);
  }

  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->SubVector(i_surf_)->Data()->Mesh();
  int ncells_surf = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(i_domain_)->Data();
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);
  Epetra_MultiVector& domain_Pu_f = *domain_Pu->ViewComponent(face_entity, false);

  // Approach 2
  double damp = 1.;
  if (damp_the_spurt_) {
    for (int cs = 0; cs != ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, cs);
      double p_old = GetDomainFaceValue(*u->SubVector(i_domain_)->Data(), f);
      double p_Pu = GetDomainFaceValue(*Pu->SubVector(i_domain_)->Data(), f);
      double p_new = p_old - p_Pu;
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        double my_damp = ((patm + cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SPURT (sc=" << surf_mesh->cell_map(false).GID(cs)
                    << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp
                    << std::endl;
      }
    }

    double proc_damp = damp;
    domain_Pu_f.Comm().MinAll(&proc_damp, &damp, 1);
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
int
MPCDelegateWater::ModifyCorrection_WaterSpurtCap(double h,
                                                 Teuchos::RCP<const TreeVector> res,
                                                 Teuchos::RCP<const TreeVector> u,
                                                 Teuchos::RCP<TreeVector> Pu,
                                                 double damp)
{
  const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->SubVector(i_surf_)->Data()->Mesh();
  int ncells_surf = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  Teuchos::RCP<const CompositeVector> domain_u = u->SubVector(i_domain_)->Data();
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();

  // Approach 3
  int n_modified = 0;
  if (cap_the_spurt_) {
    for (int cs = 0; cs != ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, cs);

      double p_old = GetDomainFaceValue(*u->SubVector(i_domain_)->Data(), f);
      double p_Pu = GetDomainFaceValue(*Pu->SubVector(i_domain_)->Data(), f);
      double p_new = p_old - p_Pu / damp;
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        double p_corrected = p_old - (patm + cap_size_);
        SetDomainFaceValue(*domain_Pu, f, p_corrected);

        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  CAPPING THE SPURT (sc=" << surf_mesh->cell_map(false).GID(cs)
                    << ",f=" << u->SubVector(i_domain_)->Data()->Mesh()->face_map(false).GID(f)
                    << "): p_old = " << p_old << ", p_new = " << p_new
                    << ", p_capped = " << p_old - p_corrected << std::endl;
      } else if ((p_new < patm) && (p_old > patm)) {
        // strange attempt to kick NKA when it goes back under?
        // double p_corrected = p_old - (patm - cap_size_);
        // SetDomainFaceValue(*domain_Pu, f, p_corrected);
        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  INVERSE SPURT (sc=" << surf_mesh->cell_map(false).GID(cs)
                    << "): p_old = " << p_old << ", p_new = " << p_new << std::endl;
      }
    }
  }
  return n_modified;
}


// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
double
MPCDelegateWater::ModifyCorrection_SaturatedSpurtDamp(double h,
                                                      Teuchos::RCP<const TreeVector> res,
                                                      Teuchos::RCP<const TreeVector> u,
                                                      Teuchos::RCP<TreeVector> Pu)
{
  const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh = u->SubVector(i_domain_)->Data()->Mesh();
  int ncells_domain = domain_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  const Epetra_MultiVector& domain_p_c =
    *u->SubVector(i_domain_)->Data()->ViewComponent("cell", false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);

  // Approach 2
  double damp = 1.;
  if (damp_the_sat_spurt_) {
    for (int c = 0; c != ncells_domain; ++c) {
      double p_old = domain_p_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c];
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        double my_damp = ((patm + cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SATURATED SPURT (c=" << domain_mesh->cell_map(false).GID(c)
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
int
MPCDelegateWater::ModifyCorrection_SaturatedSpurtCap(double h,
                                                     Teuchos::RCP<const TreeVector> res,
                                                     Teuchos::RCP<const TreeVector> u,
                                                     Teuchos::RCP<TreeVector> Pu,
                                                     double damp)
{
  const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh = u->SubVector(i_domain_)->Data()->Mesh();
  int ncells_domain = domain_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  const Epetra_MultiVector& domain_p_c =
    *u->SubVector(i_domain_)->Data()->ViewComponent("cell", false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);

  // Approach 3
  int n_modified = 0;
  if (cap_the_sat_spurt_) {
    for (int c = 0; c != ncells_domain; ++c) {
      double p_old = domain_p_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c] / damp;
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        domain_Pu_c[0][c] = p_old - (patm + cap_size_);
        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  CAPPING THE SATURATED SPURT (c=" << domain_mesh->cell_map(false).GID(c)
                    << "): p_old = " << p_old << ", p_new = " << p_new
                    << ", p_capped = " << p_old - domain_Pu_c[0][c] << std::endl;
      }
    }
  }

  return n_modified;
}

// Approach 2: damping of the spurt -- limit the max oversaturated pressure
//  using a global damping term.
double
MPCDelegateWater::ModifyCorrection_DesaturatedSpurtDamp(double h,
                                                        Teuchos::RCP<const TreeVector> res,
                                                        Teuchos::RCP<const TreeVector> u,
                                                        Teuchos::RCP<TreeVector> Pu)
{
  const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh = u->SubVector(i_domain_)->Data()->Mesh();
  int ncells_domain = domain_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  const Epetra_MultiVector& domain_p_c =
    *u->SubVector(i_domain_)->Data()->ViewComponent("cell", false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);

  // Approach 2
  double damp = 1.;
  if (damp_the_desat_spurt_) {
    for (int c = 0; c != ncells_domain; ++c) {
      double p_old = domain_p_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c];
      if ((p_new < patm - cap_size_) && (p_old >= patm)) {
        double my_damp = ((patm - cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE DESATURATED SPURT (c=" << domain_mesh->cell_map(false).GID(c)
                    << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp
                    << std::endl;
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
int
MPCDelegateWater::ModifyCorrection_DesaturatedSpurtCap(double h,
                                                       Teuchos::RCP<const TreeVector> res,
                                                       Teuchos::RCP<const TreeVector> u,
                                                       Teuchos::RCP<TreeVector> Pu,
                                                       double damp)
{
  const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh = u->SubVector(i_domain_)->Data()->Mesh();
  int ncells_domain = domain_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  const Epetra_MultiVector& domain_p_c =
    *u->SubVector(i_domain_)->Data()->ViewComponent("cell", false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->SubVector(i_domain_)->Data();
  Epetra_MultiVector& domain_Pu_c = *domain_Pu->ViewComponent("cell", false);

  // Approach 3
  int n_modified = 0;
  if (cap_the_desat_spurt_) {
    for (int c = 0; c != ncells_domain; ++c) {
      double p_old = domain_p_c[0][c];
      double p_new = p_old - domain_Pu_c[0][c] / damp;
      if ((p_new < patm - cap_size_) && (p_old >= patm)) {
        domain_Pu_c[0][c] = p_old - (patm - cap_size_);
        n_modified++;
        if (vo_->os_OK(Teuchos::VERB_HIGH))
          std::cout << "  CAPPING THE DESATURATED SPURT (c=" << domain_mesh->cell_map(false).GID(c)
                    << "): p_old = " << p_old << ", p_new = " << p_new
                    << ", p_capped = " << p_old - domain_Pu_c[0][c] << std::endl;
      }
    }
  }
  return n_modified;
}


// modify predictor via heuristic stops spurting in the surface flow
bool
MPCDelegateWater::ModifyPredictor_Heuristic(double h, const Teuchos::RCP<TreeVector>& u)
{
  bool modified = false;
  if (modify_predictor_heuristic_) {
    if (vo_->os_OK(Teuchos::VERB_HIGH))
      *vo_->os() << "  MPCWaterCoupler: Modifying predictor with water heuristic" << std::endl;

    Teuchos::RCP<CompositeVector> domain_u = u->SubVector(i_domain_)->Data();
    Epetra_MultiVector& surf_u_c = *u->SubVector(i_surf_)->Data()->ViewComponent("cell", false);

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->SubVector(i_surf_)->Data()->Mesh();

    const Epetra_MultiVector& surf_u_prev_c =
      *S_->Get<CompositeVector>("surface-pressure", tag_current_).ViewComponent("cell", false);
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);
    int ncells = surf_u_c.MyLength();
    for (int c = 0; c != ncells; ++c) {
      int f = surf_mesh->entity_get_parent(AmanziMesh::CELL, c);

      double dp = surf_u_c[0][c] - surf_u_prev_c[0][c];
      double pnew = surf_u_c[0][c] - patm;
      double pold = surf_u_prev_c[0][c] - patm;

      if (pnew > 0) {
        if (dp > pnew) {
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            *vo_->os() << "CHANGING (first over?): p = " << surf_u_c[0][c] << " to "
                       << patm + cap_size_ << std::endl;
          surf_u_c[0][c] = patm + cap_size_;
          SetDomainFaceValue(*domain_u, f, surf_u_c[0][c]);

        } else if (pold > 0 && dp > pold) {
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            *vo_->os() << "CHANGING (second over?): p = " << surf_u_c[0][c] << " to "
                       << patm + 2 * pold << std::endl;
          surf_u_c[0][c] = patm + 2 * pold;
          SetDomainFaceValue(*domain_u, f, surf_u_c[0][c]);
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
bool
MPCDelegateWater::ModifyPredictor_WaterSpurtDamp(double h, const Teuchos::RCP<TreeVector>& u)
{
  // Approach 2
  if (modify_predictor_spurt_damping_) {
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->SubVector(i_surf_)->Data()->Mesh();
    int ncells_surf = surf_mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

    Epetra_MultiVector& surf_pnew_c = *u->SubVector(i_surf_)->Data()->ViewComponent("cell", false);
    Teuchos::RCP<CompositeVector> domain_pnew = u->SubVector(i_domain_)->Data();
    Key key_ss = Keys::getKey(domain_ss_, "pressure");


    auto domain_pold = S_->GetPtr<CompositeVector>(key_ss, tag_current_);

    int rank = surf_mesh->get_comm()->MyPID();
    double damp = 1.;
    for (unsigned int cs = 0; cs != ncells_surf; ++cs) {
      AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, cs);
      double p_old = GetDomainFaceValue(*domain_pold, f);
      double p_new = GetDomainFaceValue(*domain_pnew, f);
      if ((p_new > patm + cap_size_) && (p_old < patm)) {
        // first over
        double my_damp = ((patm + cap_size_) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SPURT (1st over) (sc=" << surf_mesh->cell_map(false).GID(cs)
                    << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp
                    << std::endl;
      } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
        // second over
        double my_damp = ((patm + 2 * (p_old - patm)) - p_old) / (p_new - p_old);
        damp = std::min(damp, my_damp);
        if (vo_->os_OK(Teuchos::VERB_EXTREME))
          std::cout << "   DAMPING THE SPURT (2nd over) (sc=" << surf_mesh->cell_map(false).GID(cs)
                    << "): p_old = " << p_old << ", p_new = " << p_new << ", coef = " << my_damp
                    << std::endl;
      }
    }

    double proc_damp = damp;
    if (domain_pnew->HasComponent("face")) {
      Epetra_MultiVector& domain_pnew_f = *domain_pnew->ViewComponent("face", false);
      domain_pnew_f.Comm().MinAll(&proc_damp, &damp, 1);
    } else if (domain_pnew->HasComponent("boundary face")) {
      Epetra_MultiVector& domain_pnew_f = *domain_pnew->ViewComponent("boundary face", false);
      domain_pnew_f.Comm().MinAll(&proc_damp, &damp, 1);
    }

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
        AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, cs);
        double p_old = GetDomainFaceValue(*domain_pold, f);
        double p_new = GetDomainFaceValue(*domain_pnew, f);
        p_new = (p_new - p_old) / damp + p_old;
        if ((p_new > patm + cap_size_) && (p_old < patm)) {
          // first over
          double new_value = patm + cap_size_;
          SetDomainFaceValue(*domain_pnew, f, new_value);
          surf_pnew_c[0][cs] = patm + new_value;
          if (vo_->os_OK(Teuchos::VERB_HIGH))
            std::cout << "  CAPPING THE SPURT (1st over) (sc=" << surf_mesh->cell_map(false).GID(cs)
                      << "): p_old = " << p_old << ", p_new = " << p_new
                      << ", p_capped = " << new_value << std::endl;
        } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
          // second over
          double new_value = patm + 2 * (p_old - patm);
          SetDomainFaceValue(*domain_pnew, f, new_value);
          surf_pnew_c[0][cs] = new_value;

          if (vo_->os_OK(Teuchos::VERB_HIGH))
            std::cout << "  CAPPING THE SPURT (2nd over) (sc=" << surf_mesh->cell_map(false).GID(cs)
                      << "): p_old = " << p_old << ", p_new = " << p_new
                      << ", p_capped = " << new_value << std::endl;
        } else {
          surf_pnew_c[0][cs] = GetDomainFaceValue(*domain_pnew, f);
        }
      }
    }

    return damp < 1.0;
  }

  return false;
}


// modify predictor via heuristic stops spurting in the surface flow
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
          AmanziMesh::Entity_ID f = surf_mesh->entity_get_parent(AmanziMesh::CELL, c);
          SetDomainFaceValue(*domain_pnew, f, surf_pnew_c[0][c]);
        }
      }
    }
  }
  return modified;
}

} // namespace Amanzi
