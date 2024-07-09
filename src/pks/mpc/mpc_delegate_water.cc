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
                                   const std::string& domain,
                                   const std::string& domain_surf)
  : plist_(plist),
    S_(S),
    i_domain_(-1),
    i_surf_(-1),
    i_Tdomain_(-1),
    i_Tsurf_(-1),
    domain_ss_(domain),
    domain_surf_(domain_surf)
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
  vo_ = Teuchos::rcp(new VerboseObject(S->GetMesh(domain)->getComm(), plist->name(), *plist_));
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
    if (Pu->getSubVector(i_domain_)->getData()->hasComponent("face")) {
      face_entity = "face";
    } else if (Pu->getSubVector(i_domain_)->getData()->hasComponent("boundary_face")) {
      face_entity = "boundary_face";
    } else {
      Errors::Message message("Subsurface vector does not have face component.");
      Exceptions::amanzi_throw(message);
    }

    auto domain_Pu_f = Pu->getSubVector(i_domain_)->getData()->viewComponent(face_entity, false);
    double face_limiter(face_limiter_);

    Kokkos::parallel_reduce(
      "MPCDelegateWater::ModifyCorrection_WaterFaceLimiter",
      domain_Pu_f.extent(0),
      KOKKOS_LAMBDA(const int& f, int& l_modified) {
        if (fabs(domain_Pu_f(f, 0)) > face_limiter) {
          domain_Pu_f(f, 0) = domain_Pu_f(f, 0) > 0. ? face_limiter : -face_limiter;
          l_modified++;
        }
      },
      n_modified);
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
  if (Pu->getSubVector(i_domain_)->getData()->hasComponent("face")) {
    face_entity = "face";
  } else if (Pu->getSubVector(i_domain_)->getData()->hasComponent("boundary_face")) {
    face_entity = "boundary_face";
  } else {
    Errors::Message message("Subsurface vector does not have face component.");
    Exceptions::amanzi_throw(message);
  }

  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->getSubVector(i_surf_)->getData()->getMesh();
  int ncells_surf =
    surf_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  Teuchos::RCP<const CompositeVector> domain_u = u->getSubVector(i_domain_)->getData();
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->getSubVector(i_domain_)->getData();
  auto domain_Pu_c = domain_Pu->viewComponent("cell", false);
  auto domain_Pu_f = domain_Pu->viewComponent(face_entity, false);
  const auto domain_u_f = domain_u->viewComponent(face_entity, false);
  auto parent_ents = surf_mesh->getEntityParents<MemSpace_kind::DEVICE>(AmanziMesh::Entity_kind::CELL);
  double cap_size(cap_size_);

  // Approach 2
  double damp = 1.;
  if (damp_the_spurt_) {
    if (face_entity == "face") {
      Kokkos::parallel_reduce(
        "MPCDelegateWater::ModifyCorrection_WaterSpurtDamp",
        ncells_surf,
        KOKKOS_LAMBDA(const int& cs, double& l_damp) {
          AmanziMesh::Entity_ID f = parent_ents(cs);
          double p_old = domain_u_f(f, 0);
          double p_Pu = domain_Pu_f(f, 0);
          double p_new = p_old - p_Pu;
          if ((p_new > patm + cap_size) && (p_old < patm)) {
            double my_damp = ((patm + cap_size) - p_old) / (p_new - p_old);
            if (my_damp < l_damp) l_damp = my_damp;
          }
        },
        Kokkos::Min<double>(damp));
    } else {
      const AmanziMesh::MeshCache& sub_mesh = u->getSubVector(i_domain_)->getData()->getMesh()->getCache();
      Kokkos::parallel_reduce(
        "MPCDelegateWater::ModifyCorrection_WaterSpurtDamp",
        ncells_surf,
        KOKKOS_LAMBDA(const int& cs, double& l_damp) {
          AmanziMesh::Entity_ID bf =
            AmanziMesh::getFaceOnBoundaryBoundaryFace(sub_mesh, parent_ents(cs));
          double p_old = domain_u_f(bf, 0);
          double p_Pu = domain_Pu_f(bf, 0);
          double p_new = p_old - p_Pu;
          if ((p_new > patm + cap_size) && (p_old < patm)) {
            double my_damp = ((patm + cap_size) - p_old) / (p_new - p_old);
            if (my_damp < l_damp) l_damp = my_damp;
          }
        },
        Kokkos::Min<double>(damp));
    }

    double proc_damp = damp;
    Teuchos::reduceAll(*domain_u->getComm(), Teuchos::REDUCE_MIN, 1, &proc_damp, &damp);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SPURT!, coef = " << damp << std::endl;
      domain_Pu->scale(damp);
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

  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->getSubVector(i_surf_)->getData()->getMesh();
  int ncells_surf =
    surf_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  std::string face_entity;
  if (Pu->getSubVector(i_domain_)->getData()->hasComponent("face")) {
    face_entity = "face";
  } else if (Pu->getSubVector(i_domain_)->getData()->hasComponent("boundary_face")) {
    face_entity = "boundary_face";
  } else {
    Errors::Message message("Subsurface vector does not have face component.");
    Exceptions::amanzi_throw(message);
  }

  Teuchos::RCP<const CompositeVector> domain_u = u->getSubVector(i_domain_)->getData();
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->getSubVector(i_domain_)->getData();
  auto domain_Pu_c = domain_Pu->viewComponent("cell", false);
  auto domain_Pu_f = domain_Pu->viewComponent(face_entity, false);
  const auto domain_u_f = domain_u->viewComponent(face_entity, false);
  auto parent_ents = surf_mesh->getEntityParents<MemSpace_kind::DEVICE>(AmanziMesh::Entity_kind::CELL);
  double cap_size(cap_size_);

  // Approach 3
  int n_modified = 0;
  if (cap_the_spurt_) {
    if (face_entity == "face") {
      Kokkos::parallel_reduce(
        "MPCDelegateWater::ModifyCorrection_WaterSpurtCap",
        ncells_surf,
        KOKKOS_LAMBDA(const int& cs, int& l_modified) {
          AmanziMesh::Entity_ID f = parent_ents(cs);

          double p_old = domain_u_f(f, 0);
          double p_Pu = domain_Pu_f(f, 0);
          double p_new = p_old - p_Pu / damp;
          if ((p_new > patm + cap_size) && (p_old < patm)) {
            double p_corrected = p_old - (patm + cap_size);
            domain_Pu_f(f, 0) = p_corrected;

            l_modified++;

          } else if ((p_new < patm) && (p_old > patm)) {
            // strange attempt to kick NKA when it goes back under?
            // double p_corrected = p_old - (patm - cap_size);
            // domain_Pu_f(f,0) = p_corrected;
            l_modified++;
          }
        },
        n_modified);

    } else if (face_entity == "boundary_face") {
      const AmanziMesh::MeshCache& sub_mesh = u->getSubVector(i_domain_)->getData()->getMesh()->getCache();

      Kokkos::parallel_reduce(
        "MPCDelegateWater::ModifyCorrection_WaterSpurtCap",
        ncells_surf,
        KOKKOS_LAMBDA(const int& cs, int& l_modified) {
          AmanziMesh::Entity_ID f =
            AmanziMesh::getFaceOnBoundaryBoundaryFace(sub_mesh, parent_ents(cs));

          double p_old = domain_u_f(f, 0);
          double p_Pu = domain_Pu_f(f, 0);
          double p_new = p_old - p_Pu / damp;
          if ((p_new > patm + cap_size) && (p_old < patm)) {
            double p_corrected = p_old - (patm + cap_size);
            domain_Pu_f(f, 0) = p_corrected;

            l_modified++;

          } else if ((p_new < patm) && (p_old > patm)) {
            // strange attempt to kick NKA when it goes back under?
            // double p_corrected = p_old - (patm - cap_size);
            // domain_Pu_f(f,0) = p_corrected;
            l_modified++;
          }
        },
        n_modified);
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

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh =
    u->getSubVector(i_domain_)->getData()->getMesh();
  int ncells_domain =
    domain_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  const auto domain_p_c = u->getSubVector(i_domain_)->getData()->viewComponent("cell", false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->getSubVector(i_domain_)->getData();
  auto domain_Pu_c = domain_Pu->viewComponent("cell", false);
  double cap_size(cap_size_);

  // Approach 2
  double damp = 1.;
  if (damp_the_sat_spurt_) {
    Kokkos::parallel_reduce(
      "MPCDelegateWater::ModifyCorrection_SaturatedSpurtDamp",
      ncells_domain,
      KOKKOS_LAMBDA(const int& c, double& l_damp) {
        double p_old = domain_p_c(c, 0);
        double p_new = p_old - domain_Pu_c(c, 0);
        if ((p_new > patm + cap_size) && (p_old < patm)) {
          double my_damp = ((patm + cap_size) - p_old) / (p_new - p_old);
          if (my_damp < l_damp) l_damp = my_damp;
        }
      },
      Kokkos::Min<double>(damp));

    double proc_damp = damp;
    Teuchos::reduceAll(*domain_Pu->getComm(), Teuchos::REDUCE_MIN, 1, &proc_damp, &damp);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SATURATED SPURT!, coef = " << damp << std::endl;
      domain_Pu->scale(damp);
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

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh =
    u->getSubVector(i_domain_)->getData()->getMesh();
  int ncells_domain =
    domain_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  const auto domain_p_c = u->getSubVector(i_domain_)->getData()->viewComponent("cell", false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->getSubVector(i_domain_)->getData();
  auto domain_Pu_c = domain_Pu->viewComponent("cell", false);
  double cap_size(cap_size_);

  // Approach 3
  int n_modified = 0;
  if (cap_the_sat_spurt_) {
    Kokkos::parallel_reduce(
      "MPCDelegateWater::ModifyCorrection_SaturatedSpurtCap",
      ncells_domain,
      KOKKOS_LAMBDA(const int& c, int& l_modified) {
        double p_old = domain_p_c(c, 0);
        double p_new = p_old - domain_Pu_c(c, 0) / damp;
        if ((p_new > patm + cap_size) && (p_old < patm)) {
          domain_Pu_c(c, 0) = p_old - (patm + cap_size);
          l_modified++;
        }
      },
      n_modified);
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

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh =
    u->getSubVector(i_domain_)->getData()->getMesh();
  int ncells_domain =
    domain_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  const auto domain_p_c = u->getSubVector(i_domain_)->getData()->viewComponent("cell", false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->getSubVector(i_domain_)->getData();
  auto domain_Pu_c = domain_Pu->viewComponent("cell", false);
  double cap_size(cap_size_);

  // Approach 2
  double damp = 1.;
  if (damp_the_desat_spurt_) {
    Kokkos::parallel_reduce(
      "MPCDelegateWater::ModifyCorrection_DesaturatedSpurtDamp",
      ncells_domain,
      KOKKOS_LAMBDA(const int& c, double& l_damp) {
        double p_old = domain_p_c(c, 0);
        double p_new = p_old - domain_Pu_c(c, 0);
        if ((p_new < patm - cap_size) && (p_old >= patm)) {
          double my_damp = ((patm - cap_size) - p_old) / (p_new - p_old);
          if (my_damp < l_damp) l_damp = my_damp;
        }
      },
      Kokkos::Min<double>(damp));

    double proc_damp = damp;
    Teuchos::reduceAll(*domain_Pu->getComm(), Teuchos::REDUCE_MIN, 1, &proc_damp, &damp);
    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE DESATURATED SPURT!, coef = " << damp << std::endl;
      domain_Pu->scale(damp);
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

  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh =
    u->getSubVector(i_domain_)->getData()->getMesh();
  int ncells_domain =
    domain_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  const auto domain_p_c = u->getSubVector(i_domain_)->getData()->viewComponent("cell", false);
  Teuchos::RCP<CompositeVector> domain_Pu = Pu->getSubVector(i_domain_)->getData();
  auto domain_Pu_c = domain_Pu->viewComponent("cell", false);
  double cap_size(cap_size_);

  // Approach 3
  int n_modified = 0;
  if (cap_the_desat_spurt_) {
    Kokkos::parallel_reduce(
      "MPCDelegateWater::ModifyCorrection_DesaturatedSpurtCap",
      ncells_domain,
      KOKKOS_LAMBDA(const int& c, int& l_modified) {
        double p_old = domain_p_c(c, 0);
        double p_new = p_old - domain_Pu_c(c, 0) / damp;
        if ((p_new < patm - cap_size) && (p_old >= patm)) {
          domain_Pu_c(c, 0) = p_old - (patm - cap_size);
          l_modified++;
        }
      },
      n_modified);
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

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->getSubVector(i_surf_)->getData()->getMesh();
    auto parent_ents = surf_mesh->getEntityParents<MemSpace_kind::DEVICE>(AmanziMesh::Entity_kind::CELL);
    const double& patm = S_->Get<double>("atmospheric_pressure", Tags::DEFAULT);

    std::string face_entity;
    if (u->getSubVector(i_domain_)->getData()->hasComponent("face")) {
      face_entity = "face";
    } else if (u->getSubVector(i_domain_)->getData()->hasComponent("boundary_face")) {
      face_entity = "boundary_face";
    } else {
      Errors::Message message("Subsurface vector does not have face component.");
      Exceptions::amanzi_throw(message);
    }

    Key key_surf = Keys::getKey(domain_surf_, "pressure");
    const auto surf_u_prev_c =
      S_->Get<CompositeVector>(key_surf, tag_current_).viewComponent("cell", false);
    auto surf_u_c = u->getSubVector(i_surf_)->getData()->viewComponent("cell", false);
    Teuchos::RCP<CompositeVector> domain_u = u->getSubVector(i_domain_)->getData();
    auto domain_u_f = domain_u->viewComponent(face_entity, false);
    double cap_size(cap_size_);

    int ncells_surf = surf_u_c.extent(0);
    if (face_entity == "face") {
      Kokkos::parallel_for(
        "MPCDelegateWater::ModifyPredictor_Heuristic", ncells_surf, KOKKOS_LAMBDA(const int& c) {
          AmanziMesh::Entity_ID f = parent_ents(c);

          double dp = surf_u_c(c, 0) - surf_u_prev_c(c, 0);
          double pnew = surf_u_c(c, 0) - patm;
          double pold = surf_u_prev_c(c, 0) - patm;

          if (pnew > 0) {
            if (dp > pnew) {
              surf_u_c(c, 0) = patm + cap_size;
              domain_u_f(f, 0) = surf_u_c(c, 0);

            } else if (pold > 0 && dp > pold) {
              surf_u_c(c, 0) = patm + 2 * pold;
              domain_u_f(f, 0) = surf_u_c(c, 0);
            }
          }
        });
    } else {
      const AmanziMesh::MeshCache& sub_mesh = u->getSubVector(i_domain_)->getData()->getMesh().get()->getCache();
      Kokkos::parallel_for(
        "MPCDelegateWater::ModifyPredictor_Heuristic", ncells_surf, KOKKOS_LAMBDA(const int& c) {
          AmanziMesh::Entity_ID f =
            AmanziMesh::getFaceOnBoundaryBoundaryFace(sub_mesh, parent_ents(c));

          double dp = surf_u_c(c, 0) - surf_u_prev_c(c, 0);
          double pnew = surf_u_c(c, 0) - patm;
          double pold = surf_u_prev_c(c, 0) - patm;

          if (pnew > 0) {
            if (dp > pnew) {
              surf_u_c(c, 0) = patm + cap_size;
              domain_u_f(f, 0) = surf_u_c(c, 0);

            } else if (pold > 0 && dp > pold) {
              surf_u_c(c, 0) = patm + 2 * pold;
              domain_u_f(f, 0) = surf_u_c(c, 0);
            }
          }
        });
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

    Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->getSubVector(i_surf_)->getData()->getMesh();
    int ncells_surf =
      surf_mesh->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

    std::string face_entity;
    if (u->getSubVector(i_domain_)->getData()->hasComponent("face")) {
      face_entity = "face";
    } else if (u->getSubVector(i_domain_)->getData()->hasComponent("boundary_face")) {
      face_entity = "boundary_face";
    } else {
      Errors::Message message("Subsurface vector does not have face component.");
      Exceptions::amanzi_throw(message);
    }
    Key key_ss = Keys::getKey(domain_ss_, "pressure");
    auto domain_pold = S_->GetPtr<CompositeVector>(key_ss, tag_current_);

    Teuchos::RCP<CompositeVector> domain_pnew = u->getSubVector(i_domain_)->getData();
    double damp = 1.;

    {
      auto surf_pnew_c = u->getSubVector(i_surf_)->getData()->viewComponent("cell", false);
      auto domain_pold_f = domain_pold->viewComponent(face_entity, false);
      auto domain_pnew_f = domain_pnew->viewComponent(face_entity, false);
      auto parent_ents = surf_mesh->getEntityParents<MemSpace_kind::DEVICE>(AmanziMesh::Entity_kind::CELL);
      double cap_size(cap_size_);

      if (face_entity == "face") {
        Kokkos::parallel_reduce(
          "MPCDelegateWater::ModifyPredictor_WaterSpurtDamp",
          ncells_surf,
          KOKKOS_LAMBDA(const int& cs, double& l_damp) {
            AmanziMesh::Entity_ID f = parent_ents(cs);
            double p_old = domain_pold_f(f, 0);
            double p_new = domain_pnew_f(f, 0);
            if ((p_new > patm + cap_size) && (p_old < patm)) {
              // first over
              double my_damp = ((patm + cap_size) - p_old) / (p_new - p_old);
              l_damp = Kokkos::min(l_damp, my_damp);
            } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
              // second over
              double my_damp = ((patm + 2 * (p_old - patm)) - p_old) / (p_new - p_old);
              if (my_damp < l_damp) l_damp = my_damp;
            }
          },
          Kokkos::Min<double>(damp));
      } else {
        const AmanziMesh::MeshCache& sub_mesh = u->getSubVector(i_domain_)->getData()->getMesh()->getCache();
        Kokkos::parallel_reduce(
          "MPCDelegateWater::ModifyPredictor_WaterSpurtDamp",
          ncells_surf,
          KOKKOS_LAMBDA(const int& cs, double& l_damp) {
            AmanziMesh::Entity_ID f =
              AmanziMesh::getFaceOnBoundaryBoundaryFace(sub_mesh, parent_ents(cs));
            double p_old = domain_pold_f(f, 0);
            double p_new = domain_pnew_f(f, 0);
            if ((p_new > patm + cap_size) && (p_old < patm)) {
              // first over
              double my_damp = ((patm + cap_size) - p_old) / (p_new - p_old);
              l_damp = Kokkos::min(l_damp, my_damp);
            } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
              // second over
              double my_damp = ((patm + 2 * (p_old - patm)) - p_old) / (p_new - p_old);
              if (my_damp < l_damp) l_damp = my_damp;
            }
          },
          Kokkos::Min<double>(damp));
      }
    }

    double proc_damp = damp;
    Teuchos::reduceAll(*domain_pnew->getComm(), Teuchos::REDUCE_MIN, 1, &proc_damp, &damp);

    if (damp < 1.0) {
      if (vo_->os_OK(Teuchos::VERB_HIGH))
        *vo_->os() << "  DAMPING THE SPURT!, coef = " << damp << std::endl;

      // apply the damping
      db_->WriteVector("p_old", domain_pold.ptr());
      db_->WriteVector("p_new", domain_pnew.ptr());
      domain_pnew->update(1. - damp, *domain_pold, damp);
      db_->WriteVector("p_damped", domain_pnew.ptr());

      // undamp and cap the surface
      {
        auto surf_pnew_c = u->getSubVector(i_surf_)->getData()->viewComponent("cell", false);
        auto domain_pold_f = domain_pold->viewComponent(face_entity, false);
        auto domain_pnew_f = domain_pnew->viewComponent(face_entity, false);
        auto parent_ents = surf_mesh->getEntityParents<MemSpace_kind::DEVICE>(AmanziMesh::Entity_kind::CELL);
        double cap_size(cap_size_);

        if (face_entity == "face") {
          Kokkos::parallel_for(
            "MPCDelegateWater::ModifyPredictor_WaterSpurtDamp2",
            ncells_surf,
            KOKKOS_LAMBDA(const int& cs) {
              AmanziMesh::Entity_ID f = parent_ents(cs);
              double p_old = domain_pold_f(f, 0);
              double p_new = domain_pnew_f(f, 0);
              p_new = (p_new - p_old) / damp + p_old;
              if ((p_new > patm + cap_size) && (p_old < patm)) {
                // first over
                double new_value = patm + cap_size;
                domain_pnew_f(f, 0) = new_value;
                surf_pnew_c(cs, 0) = new_value;
              } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
                // second over
                double new_value = patm + 2 * (p_old - patm);
                domain_pnew_f(f, 0) = new_value;
                surf_pnew_c(cs, 0) = new_value;
              } else {
                surf_pnew_c(cs, 0) = domain_pnew_f(f, 0);
              }
            });
        } else {
          const AmanziMesh::MeshCache& sub_mesh = u->getSubVector(i_domain_)->getData()->getMesh()->getCache();
          Kokkos::parallel_for(
            "MPCDelegateWater::ModifyPredictor_WaterSpurtDamp2",
            ncells_surf,
            KOKKOS_LAMBDA(const int& cs) {
              AmanziMesh::Entity_ID f =
                AmanziMesh::getFaceOnBoundaryBoundaryFace(sub_mesh, parent_ents(cs));
              double p_old = domain_pold_f(f, 0);
              double p_new = domain_pnew_f(f, 0);
              p_new = (p_new - p_old) / damp + p_old;
              if ((p_new > patm + cap_size) && (p_old < patm)) {
                // first over
                double new_value = patm + cap_size;
                domain_pnew_f(f, 0) = new_value;
                surf_pnew_c(cs, 0) = new_value;
              } else if ((p_old > patm) && (p_new - p_old > p_old - patm)) {
                // second over
                double new_value = patm + 2 * (p_old - patm);
                domain_pnew_f(f, 0) = new_value;
                surf_pnew_c(cs, 0) = new_value;
              } else {
                surf_pnew_c(cs, 0) = domain_pnew_f(f, 0);
              }
            });
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
  AMANZI_ASSERT(false);
  bool modified = false;
  // if (modify_predictor_tempfromsource_) {
  //   modified = true;
  //   if (vo_->os_OK(Teuchos::VERB_HIGH))
  //     *vo_->os()
  //       << "  MPCWaterCoupler: Modifying predictor, taking surface temperature from source."
  //       << std::endl;

  //   Epetra_MultiVector& surf_Tnew_c = *u->getSubVector(i_Tsurf_)->getData()->viewComponent("cell", false);
  //   Epetra_MultiVector& surf_pnew_c = *u->getSubVector(i_surf_)->getData()->viewComponent("cell", false);


  //   // Epetra_MultiVector& domain_pnew_f = *u->getSubVector(i_domain_)->getData()
  //   //     ->viewComponent("face",false);
  //   Teuchos::RCP<CompositeVector> domain_pnew = u->getSubVector(i_domain_)->getData();

  //   const Epetra_MultiVector& Told =
  //     *S_->GetPtr<CompositeVector>("surface-temperature", tag_current_)
  //        ->viewComponent("cell", false);
  //   const Epetra_MultiVector& Tsource =
  //     *S_->GetPtr<CompositeVector>("surface-water_source_temperature", tag_next_)
  //        ->viewComponent("cell", false);
  //   const Epetra_MultiVector& hold =
  //     *S_->GetPtr<CompositeVector>("surface-ponded_depth", tag_current_)
  //        ->viewComponent("cell", false);
  //   const Epetra_MultiVector& dhsource =
  //     *S_->GetPtr<CompositeVector>("surface-water_source", tag_current_)
  //        ->viewComponent("cell", false);

  //   Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh = u->getSubVector(i_surf_)->getData()->getMesh();

  //   for (unsigned int c = 0; c != surf_Tnew_c.MyLength(); ++c) {
  //     if (surf_Tnew_c(c,0) < 271.15) {
  //       // frozen, modify predictor to ensure surface is ready to accept ice
  //       if (surf_pnew_c(c,0) < 101325.) {
  //         surf_pnew_c(c,0) = 101325.1;
  //         AmanziMesh::Entity_ID f = surf_mesh->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
  //         SetDomainFaceValue(*domain_pnew, f, surf_pnew_c(c,0));
  //       }
  //     }
  //   }
  // }
  return modified;
}

} // namespace Amanzi
