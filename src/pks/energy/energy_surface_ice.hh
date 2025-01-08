/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

//! An advection-diffusion equation for surface energy in two phases.
/*!

This is simply a surface energy equation that places a few more requirements
on the base class.  It could probably go away if we refactor to remove
hard-coded evaluators.

.. _energy_surface_ice_pk-spec:
.. admonition:: energy_surface_ice_pk-spec

    These are typically not set by the user:

    * `"coupled to subsurface via temperature`" ``[bool]`` **false** A coupling
      scheme, provided by MPC.

    * `"coupled to subsurface via flux`" ``[bool]`` **false** A coupling
      scheme, provided by MPC.

    * `"subsurface domain name`" ``[string]`` **optional** If one of the above
      coupling schemes is turned on, we need to know the subsurface mesh.
      Provided by MPC.

    INCLUDES:

    - ``[energy_pk-spec]``  See `Energy Base PK`_

*/


#ifndef PKS_ENERGY_SURFACE_ICE_HH_
#define PKS_ENERGY_SURFACE_ICE_HH_

#include "PK_Factory.hh"
#include "energy_base.hh"

namespace Amanzi {

// forward declarations
class Function;

namespace Energy {

class EnergySurfaceIce : public EnergyBase {
 public:
  EnergySurfaceIce(Teuchos::ParameterList& FElist,
                   const Teuchos::RCP<Teuchos::ParameterList>& plist,
                   const Teuchos::RCP<State>& S,
                   const Teuchos::RCP<TreeVector>& solution);

  // -- Initialize owned (dependent) variables.
  virtual void Initialize() override;

 protected:
  // -- setup the evaluators
  virtual void SetupPhysicalEvaluators_() override;

  // -- get enthalpy as a function of Dirichlet boundary data.  Note that this
  //    will get replaced by a better system when we get maps on the boundary
  //    faces.
  //  virtual void ApplyDirichletBCsToEnthalpy_(const Teuchos::Ptr<State>& S);

  virtual void AddSources_(const Tag& tag, const Teuchos::Ptr<CompositeVector>& f) override;
  virtual void AddSourcesToPrecon_(double h) override;

 protected:
  // simple heat condution term, q = K_s2a * (Tair - Tsurf)
  // air temperature function of time (not space)
  double K_surface_to_air_;

  // flags
  bool standalone_mode_;
  bool is_energy_source_term_;
  bool is_water_source_term_;
  bool is_air_conductivity_;

  // keys
  Key domain_ss_;

 private:
  // factory registration
  static RegisteredPKFactory<EnergySurfaceIce> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
