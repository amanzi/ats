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
                                   const std::string& domain_ss,
                                   const std::string& domain_surf)
  : plist_(plist),
    S_(S),
    i_domain_(-1),
    i_surf_(-1),
    i_Tdomain_(-1),
    i_Tsurf_(-1),
    domain_ss_(domain_ss),
    domain_surf_(domain_surf),
    domain_mesh_(S->GetMesh(domain_ss)),
    surf_mesh_(S->GetMesh(domain_surf))
{
  AMANZI_ASSERT(domain_ss_ != domain_surf_);
  AMANZI_ASSERT(surf_mesh_->parent() == domain_mesh_);

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

} // namespace Amanzi
