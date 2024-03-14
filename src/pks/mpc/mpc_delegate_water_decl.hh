/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Globalization for nonlinearity around the appearance/disappearance of surface water.
#ifndef AMANZI_MPC_DELEGATE_WATER_HH_
#define AMANZI_MPC_DELEGATE_WATER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "Debugger.hh"
#include "TreeVector.hh"
#include "CompositeVector.hh"
#include "State.hh"

/*!

The water delegate works to deal with discontinuities/strong nonlinearities
when surface cells shift from dry to wet (i.e. the surface pressure goes from <
atmospheric pressure to > atmospheric pressure.

These methods work to alter the predictor around this nonlinearity.

.. _mpc-delegate-water-spec:
.. admonition:: mpc-delegate-water-spec

    * `"modify predictor with heuristic`" ``[bool]`` **false** This simply
      limits the prediction to backtrack to just above atmospheric on both the
      first and second timesteps that take us over atmospheric.

    * `"modify predictor damp and cap the water spurt`" ``[bool]`` **false** The
      second both limits (caps) and damps all surface cells to ensure that all
      nearby cells are also not overshooting.  This is the preferred method.

    These methods work to alter the preconditioned correction for the same
    reasons described above.

    * `"global water face limiter`" ``[double]`` **1.e99** This is simply a limit
      to the maximum allowed size of the correction (in [Pa]) on all faces.  Any
      correction larger than this is set to this.

    * `"cap the water spurt`" ``[bool]`` **false** If a correction takes the
      pressure on a surface cell from below atmospheric (dry) to above (wet),
      the correction is set to a value which results in the new iterate to being
      CAP_SIZE over atmospheric.

    * `"damp the water spurt`" ``[bool]`` **false** A damping factor (less than
      one) is calculated to multiply the correction such that the largest
      correction takes a cell to just above atmospheric.  All faces (globally)
      are affected.

    * `"damp and cap the water spurt`" ``[bool]`` **false** None of the above
      should really be used.  Capping, when the cap is particularly severe,
      results in faces whose values are very out of equilibrium with their
      neighboring cells which are not capped.  Damping results in a tiny
      timestep in which, globally, at MOST one face can go from wet to dry.
      This looks to do a combination, in which all things are damped, but faces
      that are initially expected to go from dry to wet are pre-scaled to
      ensure that, when damped, they are also (like the biggest change) allowed
      to go from dry to wet (so that multiple cells can wet in the same step).
      This is the preferred method.

    In these methods, the following parameters are useful:

    * `"cap over atmospheric`" ``[double]`` **100** This sets the max size over
      atmospheric to which things are capped or damped. `[Pa]`

 */


namespace Amanzi {

class MPCDelegateWater {
 public:
  MPCDelegateWater(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const std::string& domain_ss,
                   const std::string& domain_surf);

  void set_db(const Teuchos::RCP<Debugger>& db) { db_ = db; }
  void set_meshes(const Teuchos::RCP<const AmanziMesh::Mesh>& domain_mesh,
                  const Teuchos::RCP<const AmanziMesh::Mesh>& surf_mesh)
  {
    domain_mesh_ = domain_mesh;
    surf_mesh_ = surf_mesh;
  }

  void set_tags(const Tag& tag_current, const Tag& tag_next)
  {
    tag_current_ = tag_current;
    tag_next_ = tag_next;
  }

  void set_indices(int i_pdomain, int i_psurf)
  {
    i_domain_ = i_pdomain;
    i_surf_ = i_psurf;
  }

  void set_indices(int i_pdomain, int i_psurf, int i_Tdomain, int i_Tsurf)
  {
    i_domain_ = i_pdomain;
    i_surf_ = i_psurf;
    i_Tdomain_ = i_Tdomain;
    i_Tsurf_ = i_Tsurf;
  }

  template <AmanziMesh::Entity_kind FaceEntity>
  bool ModifyPredictor_Heuristic(double h, const Teuchos::RCP<TreeVector>& u);

  template <AmanziMesh::Entity_kind FaceEntity>
  bool ModifyPredictor_WaterSpurtDamp(double h, const Teuchos::RCP<TreeVector>& u);

  template <AmanziMesh::Entity_kind FaceEntity>
  bool ModifyPredictor_TempFromSource(double h, const Teuchos::RCP<TreeVector>& u);

  template <AmanziMesh::Entity_kind FaceEntity>
  int ModifyCorrection_WaterFaceLimiter(double h,
                                        Teuchos::RCP<const TreeVector> res,
                                        Teuchos::RCP<const TreeVector> u,
                                        Teuchos::RCP<TreeVector> du);

  template <AmanziMesh::Entity_kind FaceEntity>
  double ModifyCorrection_WaterSpurtDamp(double h,
                                         Teuchos::RCP<const TreeVector> res,
                                         Teuchos::RCP<const TreeVector> u,
                                         Teuchos::RCP<TreeVector> du);

  template <AmanziMesh::Entity_kind FaceEntity>
  int ModifyCorrection_WaterSpurtCap(double h,
                                     Teuchos::RCP<const TreeVector> res,
                                     Teuchos::RCP<const TreeVector> u,
                                     Teuchos::RCP<TreeVector> du,
                                     double damping);

  template <AmanziMesh::Entity_kind FaceEntity>
  double ModifyCorrection_DesaturatedSpurtDamp(double h,
                                               Teuchos::RCP<const TreeVector> res,
                                               Teuchos::RCP<const TreeVector> u,
                                               Teuchos::RCP<TreeVector> du);

  template <AmanziMesh::Entity_kind FaceEntity>
  int ModifyCorrection_DesaturatedSpurtCap(double h,
                                           Teuchos::RCP<const TreeVector> res,
                                           Teuchos::RCP<const TreeVector> u,
                                           Teuchos::RCP<TreeVector> du,
                                           double damping);

  template <AmanziMesh::Entity_kind FaceEntity>
  double ModifyCorrection_SaturatedSpurtDamp(double h,
                                             Teuchos::RCP<const TreeVector> res,
                                             Teuchos::RCP<const TreeVector> u,
                                             Teuchos::RCP<TreeVector> du);

  template <AmanziMesh::Entity_kind FaceEntity>
  int ModifyCorrection_SaturatedSpurtCap(double h,
                                         Teuchos::RCP<const TreeVector> res,
                                         Teuchos::RCP<const TreeVector> u,
                                         Teuchos::RCP<TreeVector> du,
                                         double damping);

 protected:
  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<VerboseObject> vo_;
  Teuchos::RCP<Debugger> db_;

  // states
  Teuchos::RCP<State> S_;
  Tag tag_current_;
  Tag tag_next_;

  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> domain_mesh_;

  // predictor fixes
  bool modify_predictor_heuristic_;
  bool modify_predictor_spurt_damping_;
  bool modify_predictor_tempfromsource_;

  // Preconditioned correction alteration
  // -- control
  bool cap_the_spurt_;
  bool damp_the_spurt_;
  bool cap_the_sat_spurt_;
  bool damp_the_sat_spurt_;
  bool cap_the_desat_spurt_;
  bool damp_the_desat_spurt_;

  // -- parameters
  double cap_size_;
  double face_limiter_;

  // indices into the TreeVector
  int i_surf_;
  int i_domain_;
  int i_Tsurf_;
  int i_Tdomain_;

  Key domain_ss_;
  Key domain_surf_;
};

} // namespace Amanzi


#endif
