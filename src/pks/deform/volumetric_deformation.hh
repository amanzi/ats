/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Markus Berndt
           Daniil Svyatskiy
*/

//! Subsidence through bulk ice loss and cell volumetric change.
/*!

This process kernel provides for going from a cell volumetric change to an
updated unstructured mesh, and can be coupled sequentially with flow to solve
problems of flow in a subsiding porous media.

Note that this PK is slaved to the flow PK.  This PK must be advanced first,
and be weakly coupled to the flow PK (or an MPC that advances the flow PK), and
the timestep of this PK must match that of the flow PK (e.g. do not try to
subcyle one or the other).  This uses a rather hacky, unconventional use of
time tags, where we use saturations and porosities at the NEXT time, but ASSUME
they are actually the values of the CURRENT time.  This saves having to stash a
copy of these variables at the CURRENT time, which would otherwise not be used.

Note that all deformation here is vertical, and we assume that the subsurface
mesh is **perfectly columnar** and that the "build columns" parameter has been
given to the subsurface mesh.  See the Mesh_ spec for more.

The process here is governed through two options, the "deformation mode" and
the "deformation strategy."

The deformation mode describes how the cell volume change is calculated.  There
are three options here:

- "prescribed" uses a function to precribe the volume changes as a function of (t,x,y,z).

- "structural" decreases the cell volume if the porosity is above a prescribed
  "structurally connected matrix" porosity.  Think of this as bulk ice
  "propping up" the soil grains -- as that bulk ice melts, it reduces porosity
  toward the porosity in at which grains start to touch again and can be
  structurally sound.

- "saturation" is a heuristic that considers the liquid saturation directly,
  and tries to relax the liquid saturation back toward a value that is
  consistent with what the thawed soil should be.

.. todo: Move this into an evaluator!

The deformation strategy describes how the cell volume change is turned into
node coordinate changes.  Three options are available:

- "average" simply takes the average of volume change/surface area and
  horizontally averages this quantity across all neighbors.  While this has the
  advantage of being simple, it has issues when thaw gradients in the
  horizontal are not zero, as it may result in the loss of volume in a fully
  frozen cell, blowing up the pressure and breaking the code.  This is great
  when it works, but it almost never works in real problems, except in
  column-based models, where it is perfect.

- "mstk implementation" MSTK implements an iterated, local optimization method
  that, one-at-a-time, moves nodes to try and match the volumes.  This has
  fewer issues with overfitting, but doesn't always do sane things, and can be
  expensive if iterations don't work well.  This is not particularly robust
  either, but it seems to be the preferred method for 2D/3D problems.

- "global optimization" attempts to directly form and solve the minimization
  problem to find the nodal changes that result in the target volumetric
  changes.  Note this has issues with overfitting, so penalty methods are used
  to smooth the solution of the problem.  This is currently disabled.

NOTE: all deformation options are treated EXPLICITLY, and depend only upon
values from the old time.

.. _volumetric-deformation-pk-spec:
.. admonition:: volumetric-deformation-pk-spec

    * `"max timestep [s]`" ``[double]`` **inf** Sets a maximum timestep size.

    * `"deformation mode`" ``[string]`` **prescribed** See above for
      descriptions.  One of: `"prescribed`", `"structural`", or `"saturation`".

    * `"deformation strategy`" ``[string]`` **global optimization** See above
      for descriptions.  One of `"average`", `"global optimization`", or `"mstk
      implementation`"

    * `"domain name`" ``[string]`` **domain**  The mesh to deform.

    * `"surface domain name`" ``[string]`` **surface** The surface mesh.

    * `"deformation function`" ``[function-spec]`` **optional** Only used if
      "deformation mode" == "prescribed"

    EVALUATORS:

    - `"saturation_ice`" **DOMAIN-saturation_ice**
    - `"saturation_liquid`" **DOMAIN-saturation_liquid**
    - `"saturation_gas`" **DOMAIN-saturation_gas**
    - `"porosity`" **DOMAIN-porosity**
    - `"cell volume`" **DOMAIN-cell_volume**

    INCLUDES:

    - ``[pk-physical-default-spec]``

*/

#ifndef PKS_VOLUMETRIC_DEFORMATION_HH_
#define PKS_VOLUMETRIC_DEFORMATION_HH_

#include "CompositeMatrix.hh"
#include "CompositeVectorFunction.hh"
#include "Function.hh"

#include "PK.hh"
#include "PK_Factory.hh"
#include "PK_Physical.hh"
#include "PK_Helpers.hh"
//#include "MatrixVolumetricDeformation.hh"

namespace Amanzi {
namespace Deform {

class VolumetricDeformation : public PK_Physical {
 public:
  VolumetricDeformation(Teuchos::ParameterList& pk_tree,
                        const Teuchos::RCP<Teuchos::ParameterList>& global_plist,
                        const Teuchos::RCP<State>& S,
                        const Teuchos::RCP<TreeVector>& solution);

  void parseParameterList() override;

  // -- Setup data
  virtual void Setup() override;

  // -- Initialize owned (dependent) variables.
  virtual void Initialize() override;

  // -- Commit any secondary (dependent) variables.
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void FailStep(double t_old, double t_new, const Tag& tag) override;

  // -- Update diagnostics for vis.
  virtual void CalculateDiagnostics(const Tag& tag) override {}

  // -- advance via one of a few methods
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;

  virtual double get_dt() override { return dt_max_; }

  virtual void set_dt(double dt) override {}

 private:
  // strategy for calculating nodal deformation given change in cell volume
  enum DeformStrategy {
    DEFORM_STRATEGY_GLOBAL_OPTIMIZATION,
    DEFORM_STRATEGY_MSTK,
    DEFORM_STRATEGY_AVERAGE
  };
  DeformStrategy strategy_;

  // strategy for calculating change in cell volume
  enum DeformMode { DEFORM_MODE_DVDT, DEFORM_MODE_SATURATION, DEFORM_MODE_STRUCTURAL };
  DeformMode deform_mode_;
  double overpressured_limit_;

  // region that may deform
  std::string deform_region_;

  // parameters for DEFORM_MODE_DVDT
  Teuchos::RCP<Functions::CompositeVectorFunction> deform_func_;

  // parameters for DEFORM_MODE_SATURATION
  double min_porosity_;
  double deform_scaling_;
  double min_S_liq_;

  // parameters for DEFORM_MODE_STRUCTURAL
  double time_scale_, structural_vol_frac_;

  // PK timestep control
  double dt_max_;

  // meshes
  Key domain_surf_, domain_surf_3d_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf_mesh_;
  Teuchos::RCP<const AmanziMesh::Mesh> surf3d_mesh_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_nc_;
  Teuchos::RCP<AmanziMesh::Mesh> surf_mesh_nc_;
  Teuchos::RCP<AmanziMesh::Mesh> surf3d_mesh_nc_;

  // keys
  Key sat_liq_key_, sat_gas_key_, sat_ice_key_;
  Key cv_key_, del_cv_key_;
  Key poro_key_;
  Key vertex_loc_key_, vertex_loc_surf3d_key_;
  Key nodal_dz_key_, face_above_dz_key_;
  bool deformed_this_step_;

  // factory registration
  static RegisteredPKFactory<VolumetricDeformation> reg_;
};

} // namespace Deform
} // namespace Amanzi

#endif
